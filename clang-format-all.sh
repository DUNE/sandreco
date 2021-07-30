<<<<<<< HEAD
#!/bin/bash
find ./src ./include -type f -not -name '*.cmd' -and -name '*.c*' -or -name '*.h*' | xargs clang-format -i --style=file

=======
# clang-format
find ./src ./include -type f -not -name '*.cmd' -and -name '*.c*' -or -name '*.h*' | xargs clang-format -i --style=file

# refresh line endings
git add --renormalize .

>>>>>>> 97a935740b279c0d02e757f514dd1a9d67cd7dd2
