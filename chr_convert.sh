#chr_convert.sh
#By Beth Rowan
#14.09.15

#subsitutes CHRX with the chrnumber

#!/bin/sh

file=$1

sed '
s/Chr1/1/g
s/Chr2/2/g
s/Chr3/3/g
s/Chr4/4/g
s/Chr5/5/g'< ${1} > input.out

