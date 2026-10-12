//! GPU courtesy rules: when a running simulation must leave the GPU alone.
//!
//! Two kinds of rule, checked before every 64th step:
//!
//! - Pause file: while any listed file exists the run sleeps in place and
//!   resumes when it is gone. Files: `cfm_pause` in the project root, or the
//!   `;`-separated paths in CFM_PAUSE_FILE.
//! - Guards: a weekly time window and a list of programs that must have the
//!   machine to themselves. While a guard is active the run either gives up
//!   (exit code `YIELD_EXIT`, freeing the GPU memory; a runner starts it
//!   again later) or sleeps in place (CFM_GUARD_ACTION=pause).
//!
//! The guards are built in for the borrowed machine ISABELLE (school hours
//! 07:00-17:00 Monday to Friday, and Genshin Impact running) and cannot be
//! switched off there. Any machine can add more: CFM_GUARD_HOURS =
//! "1-5,07:00-17:00" (days 0 = Sunday .. 6, local time), CFM_GUARD_PROCS =
//! "a.exe;b.exe".

/// Exit code of a run that gave way to a guard (EX_TEMPFAIL).
pub const YIELD_EXIT: i32 = 75;

#[derive(Clone, Debug, PartialEq)]
pub struct Hours {
    pub day_from: u32,
    pub day_to: u32,
    pub min_from: u32,
    pub min_to: u32,
}

#[derive(Clone, Debug, Default)]
pub struct Rules {
    pub hours: Vec<Hours>,
    pub procs: Vec<String>,
}

/// "1-5,07:00-17:00" -> days 1..=5, minutes 420..1020.
pub fn parse_hours(s: &str) -> Option<Hours> {
    let (days, times) = s.split_once(',')?;
    let (d0, d1) = days.trim().split_once('-')?;
    let (t0, t1) = times.trim().split_once('-')?;
    let hm = |t: &str| -> Option<u32> {
        let (h, m) = t.trim().split_once(':')?;
        Some(h.parse::<u32>().ok()? * 60 + m.parse::<u32>().ok()?)
    };
    Some(Hours { day_from: d0.trim().parse().ok()?, day_to: d1.trim().parse().ok()?, min_from: hm(t0)?, min_to: hm(t1)? })
}

/// The rules in force on `host` with the given extra settings.
pub fn rules_for(host: &str, extra_hours: &str, extra_procs: &str) -> Rules {
    let mut r = Rules::default();
    if host.eq_ignore_ascii_case("ISABELLE") {
        r.hours.push(Hours { day_from: 1, day_to: 5, min_from: 7 * 60, min_to: 17 * 60 });
        r.procs.push("GenshinImpact.exe".to_string());
        r.procs.push("YuanShen.exe".to_string());
    }
    r.hours.extend(extra_hours.split(';').filter_map(parse_hours));
    r.procs.extend(extra_procs.split(';').map(|p| p.trim().to_string()).filter(|p| !p.is_empty()));
    r
}

/// Is (weekday 0 = Sunday, minute of the day) inside a window?
pub fn in_hours(h: &[Hours], weekday: u32, minute: u32) -> bool {
    h.iter().any(|w| weekday >= w.day_from && weekday <= w.day_to && minute >= w.min_from && minute < w.min_to)
}

/// Minutes from (weekday, minute) until the next guarded window opens (0
/// when inside one); None when there are no windows.
pub fn minutes_until_hours(h: &[Hours], weekday: u32, minute: u32) -> Option<u32> {
    if h.is_empty() {
        return None;
    }
    (0..7 * 24 * 60).find(|&d| {
        let t = minute + d;
        in_hours(h, (weekday + t / (24 * 60)) % 7, t % (24 * 60))
    })
}

/// Which of `procs` appears in a `tasklist /FO CSV /NH` listing.
pub fn proc_running<'a>(procs: &'a [String], tasklist_csv: &str) -> Option<&'a str> {
    let low = tasklist_csv.to_ascii_lowercase();
    procs.iter().find(|p| low.contains(&format!("\"{}\"", p.to_ascii_lowercase()))).map(|p| p.as_str())
}

#[cfg(windows)]
fn local_weekday_minute() -> Option<(u32, u32)> {
    #[repr(C)]
    #[derive(Default)]
    struct SystemTime {
        year: u16,
        month: u16,
        day_of_week: u16,
        day: u16,
        hour: u16,
        minute: u16,
        second: u16,
        millis: u16,
    }
    #[link(name = "kernel32")]
    extern "system" {
        fn GetLocalTime(t: *mut SystemTime);
    }
    let mut t = SystemTime::default();
    // SAFETY: GetLocalTime only writes the SYSTEMTIME it is handed.
    unsafe { GetLocalTime(&mut t) };
    Some((t.day_of_week as u32, t.hour as u32 * 60 + t.minute as u32))
}

#[cfg(not(windows))]
fn local_weekday_minute() -> Option<(u32, u32)> {
    None
}

fn tasklist() -> String {
    std::process::Command::new("tasklist")
        .args(["/FO", "CSV", "/NH"])
        .output()
        .map(|o| String::from_utf8_lossy(&o.stdout).into_owned())
        .unwrap_or_default()
}

fn pause_files() -> Vec<std::path::PathBuf> {
    match std::env::var("CFM_PAUSE_FILE") {
        Ok(s) if !s.trim().is_empty() => s.split(';').map(|p| std::path::PathBuf::from(p.trim())).collect(),
        _ => vec![std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("..").join("cfm_pause")],
    }
}

fn rules_here() -> Rules {
    let env = |k: &str| std::env::var(k).unwrap_or_default();
    rules_for(&env("COMPUTERNAME"), &env("CFM_GUARD_HOURS"), &env("CFM_GUARD_PROCS"))
}

/// Why the GPU must be left alone right now, if it must: a pause file, the
/// time window, or a guarded program.
pub fn blocked() -> Option<String> {
    if let Some(p) = pause_files().into_iter().find(|p| p.exists()) {
        return Some(format!("pause file {}", p.display()));
    }
    let rules = rules_here();
    if let Some(why) = guard_active(&rules) {
        return Some(why);
    }
    // CFM_GUARD_NEED_MIN: a run of that many minutes must end before the
    // next guarded window opens.
    let need: u32 = std::env::var("CFM_GUARD_NEED_MIN").ok().and_then(|q| q.trim().parse().ok()).unwrap_or(0);
    if need > 0 {
        if let Some((wd, min)) = local_weekday_minute() {
            if let Some(left) = minutes_until_hours(&rules.hours, wd, min) {
                if left < need {
                    return Some(format!("guarded hours start in {left} min, the run needs {need}"));
                }
            }
        }
    }
    None
}

fn guard_active(r: &Rules) -> Option<String> {
    if let Some((wd, min)) = local_weekday_minute() {
        if in_hours(&r.hours, wd, min) {
            return Some(format!("guarded hours (weekday {wd}, {:02}:{:02})", min / 60, min % 60));
        }
    }
    if !r.procs.is_empty() {
        if let Some(p) = proc_running(&r.procs, &tasklist()) {
            return Some(format!("{p} is running"));
        }
    }
    None
}

/// Called before each step; acts on every 64th call.
pub fn wait_while_paused() {
    use std::sync::atomic::{AtomicU32, Ordering};
    static COUNT: AtomicU32 = AtomicU32::new(0);
    if COUNT.fetch_add(1, Ordering::Relaxed) % 64 != 0 {
        return;
    }
    let files = pause_files();
    let rules = rules_here();
    let yield_now = std::env::var("CFM_GUARD_ACTION").map(|a| a != "pause").unwrap_or(true);
    let mut slept = false;
    loop {
        let reason = if let Some(p) = files.iter().find(|p| p.exists()) {
            format!("pause file {}", p.display())
        } else if let Some(why) = guard_active(&rules) {
            if yield_now {
                eprintln!("yielding the GPU: {why}");
                println!("yielded: {why}");
                std::process::exit(YIELD_EXIT);
            }
            why
        } else {
            break;
        };
        if !slept {
            eprintln!("paused ({reason})");
            slept = true;
        }
        std::thread::sleep(std::time::Duration::from_secs(if rules.procs.is_empty() { 2 } else { 10 }));
    }
    if slept {
        eprintln!("resumed");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn isabelle_rules_are_built_in() {
        let r = rules_for("isabelle", "", "");
        assert_eq!(r.hours, vec![Hours { day_from: 1, day_to: 5, min_from: 420, min_to: 1020 }]);
        assert!(r.procs.iter().any(|p| p == "GenshinImpact.exe"));
        let other = rules_for("ESTELLE", "", "");
        assert!(other.hours.is_empty() && other.procs.is_empty());
    }

    #[test]
    fn school_hours_window() {
        let r = rules_for("ISABELLE", "", "");
        // Monday 06:59 free, 07:00 guarded, 16:59 guarded, 17:00 free.
        assert!(!in_hours(&r.hours, 1, 6 * 60 + 59));
        assert!(in_hours(&r.hours, 1, 7 * 60));
        assert!(in_hours(&r.hours, 5, 16 * 60 + 59));
        assert!(!in_hours(&r.hours, 5, 17 * 60));
        // Saturday and Sunday are free all day.
        assert!(!in_hours(&r.hours, 6, 12 * 60));
        assert!(!in_hours(&r.hours, 0, 12 * 60));
    }

    #[test]
    fn time_left_before_the_window() {
        let r = rules_for("ISABELLE", "", "");
        // Sunday 23:00 -> Monday 07:00 is 8 h; Friday 17:00 -> Monday 07:00 is 62 h.
        assert_eq!(minutes_until_hours(&r.hours, 0, 23 * 60), Some(8 * 60));
        assert_eq!(minutes_until_hours(&r.hours, 5, 17 * 60), Some(62 * 60));
        assert_eq!(minutes_until_hours(&r.hours, 3, 12 * 60), Some(0));
        assert_eq!(minutes_until_hours(&[], 3, 12 * 60), None);
    }

    #[test]
    fn extra_rules_parse() {
        let r = rules_for("ESTELLE", "0-6,22:30-23:00", "notepad.exe; Game.exe");
        assert_eq!(r.hours, vec![Hours { day_from: 0, day_to: 6, min_from: 1350, min_to: 1380 }]);
        assert_eq!(r.procs, vec!["notepad.exe".to_string(), "Game.exe".to_string()]);
        assert!(parse_hours("nonsense").is_none());
    }

    #[test]
    fn process_match_is_exact_and_caseless() {
        let procs = vec!["GenshinImpact.exe".to_string(), "YuanShen.exe".to_string()];
        let list = "\"System\",\"4\",\"Services\",\"0\",\"144 K\"\r\n\"genshinimpact.exe\",\"9912\",\"Console\",\"1\",\"3,201,144 K\"\r\n";
        assert_eq!(proc_running(&procs, list), Some("GenshinImpact.exe"));
        let other = "\"NotGenshinImpact.exe\",\"1\",\"Console\",\"1\",\"1 K\"\r\n\"HYP.exe\",\"2\",\"Console\",\"1\",\"1 K\"\r\n";
        assert_eq!(proc_running(&procs, other), None);
    }

    /// Runner helper (no GPU): exits 0 when the GPU may be used now (and,
    /// with CFM_GUARD_NEED_MIN, for that many minutes), or YIELD_EXIT with
    /// the reason. `<exe> guard_check --ignored --nocapture`.
    #[test]
    #[ignore]
    fn guard_check() {
        match blocked() {
            Some(why) => {
                println!("blocked: {why}");
                std::process::exit(YIELD_EXIT);
            }
            None => println!("clear"),
        }
    }
}
