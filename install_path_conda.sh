#!/usr/bin/env bash

# Install AbRaCD PATH hooks in the currently active Conda environment.
# Run with: bash install_path_conda.sh
set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

if [[ -z "${CONDA_PREFIX:-}" || ! -d "$CONDA_PREFIX" ]]; then
    printf 'No active Conda environment detected.\n' >&2
    printf 'Activate the environment where you want AbRaCD available, then run this installer again.\n' >&2
    exit 1
fi

activate_dir="$CONDA_PREFIX/etc/conda/activate.d"
deactivate_dir="$CONDA_PREFIX/etc/conda/deactivate.d"
mkdir -p "$activate_dir" "$deactivate_dir"

cat > "$activate_dir/abracd_path.sh" <<EOF
# Added by AbRaCD install_path_conda.sh; scoped to this Conda environment.
export _ABRACD_OLD_PATH="\$PATH"
_abracd_root="$repo_root"
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
EOF

cat > "$deactivate_dir/abracd_path.sh" <<'EOF'
# Added by AbRaCD install_path_conda.sh; restore the PATH from before activation.
if [[ -n "${_ABRACD_OLD_PATH:-}" ]]; then
    export PATH="$_ABRACD_OLD_PATH"
    unset _ABRACD_OLD_PATH
fi
EOF

printf 'Installed AbRaCD PATH hooks in Conda environment: %s\n' "$CONDA_PREFIX"
printf 'Repository: %s\n' "$repo_root"
printf 'Reactivate this environment to apply the change. The PATH entries will be removed on deactivation.\n'
