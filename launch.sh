cmake CMakeLists.txt
if [ $? -eq 0 ]; then
    make
    if [ $? -ne 0 ]; then
         echo "MAKE FAILED, STOPPING" && exit 1
    fi
fi
for var in "$@"
do
    ./$var
done
