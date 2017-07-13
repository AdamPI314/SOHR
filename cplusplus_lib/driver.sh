#!/bin/sh

make run

# 1-->2
grep -rl \"iterationNumber\":\ 1, input | xargs sed -i 's/\"iterationNumber\":\ 1,/\"iterationNumber\":\ 2,/g'
make run

# 2--3>
grep -rl \"iterationNumber\":\ 2, input | xargs sed -i 's/\"iterationNumber\":\ 2,/\"iterationNumber\":\ 3,/g'
make run

# 3--4>
grep -rl \"iterationNumber\":\ 3, input | xargs sed -i 's/\"iterationNumber\":\ 3,/\"iterationNumber\":\ 4,/g'
make run

# 4--5>
grep -rl \"iterationNumber\":\ 4, input | xargs sed -i 's/\"iterationNumber\":\ 4,/\"iterationNumber\":\ 5,/g'
make run

# 5--6>
grep -rl \"iterationNumber\":\ 5, input | xargs sed -i 's/\"iterationNumber\":\ 5,/\"iterationNumber\":\ 6,/g'
make run

# 6--7>
grep -rl \"iterationNumber\":\ 6, input | xargs sed -i 's/\"iterationNumber\":\ 6,/\"iterationNumber\":\ 7,/g'
make run

# 7--8>
grep -rl \"iterationNumber\":\ 7, input | xargs sed -i 's/\"iterationNumber\":\ 7,/\"iterationNumber\":\ 8,/g'
make run

# 8--9>
grep -rl \"iterationNumber\":\ 8, input | xargs sed -i 's/\"iterationNumber\":\ 8,/\"iterationNumber\":\ 9,/g'
make run

# 9--10>
grep -rl \"iterationNumber\":\ 9, input | xargs sed -i 's/\"iterationNumber\":\ 9,/\"iterationNumber\":\ 10,/g'
make run
