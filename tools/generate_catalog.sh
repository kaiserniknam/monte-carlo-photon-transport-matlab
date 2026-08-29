#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
output="$repo_root/docs/CODE_CATALOG.md"
tmp="$(mktemp)"

{
    echo "# Code Catalog"
    echo
    echo "This catalog is generated from the MATLAB help headers. Files are clustered by scientific purpose while retaining their original chronological identifiers."
    echo
    echo "| Cluster | File | Role | Summary |"
    echo "|---|---|---|---|"
    while IFS= read -r file; do
        relative="${file#"$repo_root"/}"
        cluster="$(basename "$(dirname "$file")")"
        name="$(basename "$file")"
        role="$(sed -n '3s/^% Version role: //p' "$file")"
        summary="$(awk '
            NR <= 20 && /^%/ {
                line=$0
                sub(/^%+[[:space:]]*/, "", line)
                if (line ~ /^(Repository group|Version role|Research note|Review configuration)/ || line == "") next
                print line
                exit
            }
        ' "$file")"
        summary="${summary//$'\r'/}"
        summary="${summary//|/\\|}"
        printf '| `%s` | [`%s`](../%s) | %s | %s |\n' "$cluster" "$name" "$relative" "$role" "$summary"
    done < <(find "$repo_root/code" -type f -name '*.m' | sort -V)
} > "$tmp"

mv "$tmp" "$output"
