export XDG_RUNTIME_DIR=/tmp/adit-test2-runtime-dir
if ! test -d "${XDG_RUNTIME_DIR}"; then
    mkdir "${XDG_RUNTIME_DIR}"
    chmod 0700 "${XDG_RUNTIME_DIR}"
fi
weston-launch -- --backend=fbdev-backend.so > error.txt 2>&1
