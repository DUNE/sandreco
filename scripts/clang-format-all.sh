find ./src ./include -type f -not -name '*.cmd' -and -name '*.c*' -or -name '*.h*' | xargs clang-format -i --style=file
