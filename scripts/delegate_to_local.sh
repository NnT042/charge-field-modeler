#!/bin/bash
# delegate_to_local.sh — Send a coding task to LM Studio's local server
# Usage: ./delegate_to_local.sh "Write an impl block for SpinLevel that..."
# Or:    echo "prompt text" | ./delegate_to_local.sh -
#
# Requires: LM Studio running with server enabled (system tray → Start Server)
# Default endpoint: http://localhost:1234/v1/chat/completions

LMSTUDIO_URL="${LMSTUDIO_URL:-http://localhost:1234/v1/chat/completions}"
MAX_TOKENS="${MAX_TOKENS:-4096}"
TEMPERATURE="${TEMPERATURE:-0.2}"

# Read prompt from argument or stdin
if [ "$1" = "-" ]; then
    PROMPT=$(cat)
elif [ -n "$1" ]; then
    PROMPT="$1"
else
    echo "Usage: delegate_to_local.sh \"prompt\" or echo \"prompt\" | delegate_to_local.sh -" >&2
    exit 1
fi

# System prompt tuned for CFM code generation
SYSTEM_PROMPT="You are a code generation assistant for the Charge Field Modeler project. The stack is Godot 4 + Rust (via gdext GDExtension) + Vulkan compute shaders (GLSL). Generate clean, well-commented code. Do not explain unless asked — just output the code. Use idiomatic Rust with proper error handling. For GDScript, use static typing. For GLSL, target Vulkan 1.1 / SPIR-V."

# Escape for JSON
SYSTEM_JSON=$(printf '%s' "$SYSTEM_PROMPT" | python3 -c 'import json,sys; print(json.dumps(sys.stdin.read()))')
PROMPT_JSON=$(printf '%s' "$PROMPT" | python3 -c 'import json,sys; print(json.dumps(sys.stdin.read()))')

RESPONSE=$(curl -s -X POST "$LMSTUDIO_URL" \
    -H "Content-Type: application/json" \
    -d "{
        \"messages\": [
            {\"role\": \"system\", \"content\": $SYSTEM_JSON},
            {\"role\": \"user\", \"content\": $PROMPT_JSON}
        ],
        \"max_tokens\": $MAX_TOKENS,
        \"temperature\": $TEMPERATURE,
        \"stream\": false
    }" 2>/dev/null)

# Check if curl succeeded
if [ $? -ne 0 ]; then
    echo "ERROR: Could not connect to LM Studio at $LMSTUDIO_URL" >&2
    echo "Make sure LM Studio is running and the server is started (system tray → Start Server)" >&2
    exit 1
fi

# Extract the response text
OUTPUT=$(echo "$RESPONSE" | python3 -c '
import json, sys
try:
    data = json.load(sys.stdin)
    if "choices" in data and len(data["choices"]) > 0:
        print(data["choices"][0]["message"]["content"])
    elif "error" in data:
        print("LM Studio error: " + json.dumps(data["error"]), file=sys.stderr)
        sys.exit(1)
    else:
        print("Unexpected response format: " + json.dumps(data)[:200], file=sys.stderr)
        sys.exit(1)
except json.JSONDecodeError:
    print("ERROR: Invalid JSON response from LM Studio", file=sys.stderr)
    sys.exit(1)
')

if [ $? -ne 0 ]; then
    exit 1
fi

echo "$OUTPUT"
