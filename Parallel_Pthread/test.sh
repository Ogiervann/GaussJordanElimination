
set -e

: '
echo 2:
./a.out 2 1 1 2 1
./a.out 2 2 1 2 1
./a.out 2 1 1 2 2
./a.out 2 2 1 2 2
./a.out 2 1 1 2 3
./a.out 2 2 1 2 3
./a.out 2 1 2 2 1
./a.out 2 2 2 2 1
./a.out 2 1 2 2 2
./a.out 2 2 2 2 2
./a.out 2 1 2 2 3
./a.out 2 2 2 2 3

echo 3:

./a.out 3 1 1 3 1
./a.out 3 2 1 3 1
./a.out 3 3 1 3 1
./a.out 3 1 1 3 2
./a.out 3 2 1 3 2
./a.out 3 3 1 3 2
./a.out 3 1 1 3 3
./a.out 3 2 1 3 3
./a.out 3 3 1 3 3
./a.out 3 1 2 3 1
./a.out 3 2 2 3 1
./a.out 3 3 2 3 1
./a.out 3 1 2 3 2
./a.out 3 2 2 3 2
./a.out 3 3 2 3 2
./a.out 3 1 2 3 3
./a.out 3 2 2 3 3
./a.out 3 3 2 3 3
./a.out 3 1 3 3 1
./a.out 3 2 3 3 1
./a.out 3 3 3 3 1
./a.out 3 1 3 3 2
./a.out 3 2 3 3 2
./a.out 3 3 3 3 2
./a.out 3 1 3 3 3
./a.out 3 2 3 3 3
./a.out 3 3 3 3 3
'
echo 4:

./a.out 4 1 1 10 1
./a.out 4 2 1 10 1
./a.out 4 3 1 10 1
./a.out 4 4 1 10 1
./a.out 4 1 1 10 2
./a.out 4 2 1 10 2
./a.out 4 3 1 10 2
./a.out 4 4 1 10 2
./a.out 4 1 1 10 3
./a.out 4 2 1 10 3
./a.out 4 3 1 10 3
./a.out 4 4 1 10 3


./a.out 4 1 2 10 1
./a.out 4 2 2 10 1
./a.out 4 3 2 10 1
./a.out 4 4 2 10 1
./a.out 4 1 2 10 2
./a.out 4 2 2 10 2
./a.out 4 3 2 10 2
./a.out 4 4 2 10 2
./a.out 4 1 2 10 3
./a.out 4 2 2 10 3
./a.out 4 3 2 10 3
./a.out 4 4 2 10 3

./a.out 4 1 3 10 1
./a.out 4 2 3 10 1
./a.out 4 3 3 10 1
./a.out 4 4 3 10 1
./a.out 4 1 3 10 2
./a.out 4 2 3 10 2
./a.out 4 3 3 10 2
./a.out 4 4 3 10 2
./a.out 4 1 3 10 3
./a.out 4 2 3 10 3
./a.out 4 3 3 10 3
./a.out 4 4 3 10 3

./a.out 4 1 4 10 1
./a.out 4 2 4 10 1
./a.out 4 3 4 10 1
./a.out 4 4 4 10 1
./a.out 4 1 4 10 2
./a.out 4 2 4 10 2
./a.out 4 3 4 10 2
./a.out 4 4 4 10 2
./a.out 4 1 4 10 3
./a.out 4 2 4 10 3
./a.out 4 3 4 10 3
./a.out 4 4 4 10 3

: '
./a.out 5 1 10 1
./a.out 5 2 10 1
./a.out 5 3 10 1
./a.out 5 4 10 1
./a.out 5 5 10 1
./a.out 5 1 10 2
./a.out 5 2 10 2
./a.out 5 3 10 2
./a.out 5 4 10 2
./a.out 5 5 10 2
./a.out 5 1 10 3
./a.out 5 2 10 3
./a.out 5 3 10 3
./a.out 5 4 10 3
./a.out 5 5 10 3


./a.out 6 1 10 1
./a.out 6 2 10 1
./a.out 6 3 10 1
./a.out 6 4 10 1
./a.out 6 5 10 1
./a.out 6 6 10 1
./a.out 6 1 10 2
./a.out 6 2 10 2
./a.out 6 3 10 2
./a.out 6 4 10 2
./a.out 6 5 10 2
./a.out 6 6 10 2
./a.out 6 1 10 3
./a.out 6 2 10 3
./a.out 6 3 10 3
./a.out 6 4 10 3
./a.out 6 5 10 3
./a.out 6 6 10 3


./a.out 7 1 10 1
./a.out 7 2 10 1
./a.out 7 3 10 1
./a.out 7 4 10 1
./a.out 7 5 10 1
./a.out 7 6 10 1
./a.out 7 7 10 1
./a.out 7 1 10 2
./a.out 7 2 10 2
./a.out 7 3 10 2
./a.out 7 4 10 2
./a.out 7 5 10 2
./a.out 7 6 10 2
./a.out 7 7 10 2
./a.out 7 1 10 3
./a.out 7 2 10 3
./a.out 7 3 10 3
./a.out 7 4 10 3
./a.out 7 5 10 3
./a.out 7 6 10 3
./a.out 7 7 10 3

./a.out 8 1 10 1
./a.out 8 2 10 1
./a.out 8 3 10 1
./a.out 8 4 10 1
./a.out 8 5 10 1
./a.out 8 6 10 1
./a.out 8 7 10 1
./a.out 8 8 10 1
./a.out 8 1 10 2
./a.out 8 2 10 2
./a.out 8 3 10 2
./a.out 8 4 10 2
./a.out 8 5 10 2
./a.out 8 6 10 2
./a.out 8 7 10 2
./a.out 8 8 10 2
./a.out 8 1 10 3
./a.out 8 2 10 3
./a.out 8 3 10 3
./a.out 8 4 10 3
./a.out 8 5 10 3
./a.out 8 6 10 3
./a.out 8 7 10 3
./a.out 8 8 10 3


./a.out 9 1 10 1
./a.out 9 2 10 1
./a.out 9 3 10 1
./a.out 9 4 10 1
./a.out 9 5 10 1
./a.out 9 6 10 1
./a.out 9 7 10 1
./a.out 9 8 10 1
./a.out 9 9 10 1
./a.out 9 1 10 2
./a.out 9 2 10 2
./a.out 9 3 10 2
./a.out 9 4 10 2
./a.out 9 5 10 2
./a.out 9 6 10 2
./a.out 9 7 10 2
./a.out 9 8 10 2
./a.out 9 9 10 2
./a.out 9 1 10 3
./a.out 9 2 10 3
./a.out 9 3 10 3
./a.out 9 4 10 3
./a.out 9 5 10 3
./a.out 9 6 10 3
./a.out 9 7 10 3
./a.out 9 8 10 3
./a.out 9 9 10 3


./a.out 10 1 10 1
./a.out 10 2 10 1
./a.out 10 3 10 1
./a.out 10 4 10 1
./a.out 10 5 10 1
./a.out 10 6 10 1
./a.out 10 7 10 1
./a.out 10 8 10 1
./a.out 10 9 10 1
./a.out 10 10 10 1
./a.out 10 1 10 2
./a.out 10 2 10 2
./a.out 10 3 10 2
./a.out 10 4 10 2
./a.out 10 5 10 2
./a.out 10 6 10 2
./a.out 10 7 10 2
./a.out 10 8 10 2
./a.out 10 9 10 2
./a.out 10 10 10 2
./a.out 10 1 10 3
./a.out 10 2 10 3
./a.out 10 3 10 3
./a.out 10 4 10 3
./a.out 10 5 10 3
./a.out 10 6 10 3
./a.out 10 7 10 3
./a.out 10 8 10 3
./a.out 10 9 10 3
./a.out 10 10 10 3
'
