#!/usr/bin/env bash

# Add AbRaCD command and module directories to this user's Bash PATH.
# Run with: bash install_path_bash.sh
set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
shell_config="${HOME}/.bashrc"
start_marker="# >>> AbRaCD PATH setup >>>"
end_marker="# <<< AbRaCD PATH setup <<<"

touch "$shell_config"
temp_config="$(mktemp)"
trap 'rm -f "$temp_config"' EXIT

# Replace an existing installer-managed block, if present, to keep this
# installer repeatable and update the path if the repository has moved.
awk -v start="$start_marker" -v end="$end_marker" '
    $0 == start { removing = 1; next }
    $0 == end { removing = 0; next }
    !removing { print }
' "$shell_config" > "$temp_config"

printf -v quoted_root '%q' "$repo_root"
cat >> "$temp_config" <<EOF
$start_marker
_abracd_root=$quoted_root
for _abracd_dir in \
    "\$_abracd_root/scripts/bash" \
    "\$_abracd_root/scripts/perl" \
    "\$_abracd_root/scripts/python" \
    "\$_abracd_root/src/perl" \
    "\$_abracd_root/src/python"
do
    case ":\$PATH:" in
        *":\$_abracd_dir:"*) ;;
        *) export PATH="\$_abracd_dir:\$PATH" ;;
    esac
done
unset _abracd_root _abracd_dir
$end_marker
EOF

cp "$temp_config" "$shell_config"

printf 'Added AbRaCD script directories to %s\n' "$shell_config"
printf 'Repository: %s\n' "$repo_root"
printf 'Open a new terminal or run: source %s\n' "$shell_config"
