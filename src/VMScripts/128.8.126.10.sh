#!/bin/sh
cd /
cd Users/amos/
cd amos/

./bootstrap > /Users/amos/$1.log 2>&1
if [ $? -ne 0 ]
then
cp /Users/amos/$1.log /Users/amos/$1_Failed.log
echo "FAILED: ./bootstrap" >> /Users/amos/$1_Failed.log
/usr/bin/expect <<EOD
spawn scp /Users/amos/$1_Failed.log ssh@sauron.cs.umd.edu:VMlogs
expect "ssh@sauron.cs.umd.edu's password:"
send "123\r"
expect eof
EOD
exit
fi

./configure --with-qmake-qt4=/sw/lib/qt4-x11/bin/qmake --prefix=/usr/local/AMOS >> /Users/amos/$1.log 2>&1
if [ $? -ne 0 ]
then
cp /Users/amos/$1.log /Users/amos/$1_Failed.log
echo "FAILED: ./configure" >> /Users/amos/$1_Failed.log
/usr/bin/expect <<EOD
spawn scp /Users/amos/$1_Failed.log ssh@sauron.cs.umd.edu:VMlogs
expect "ssh@sauron.cs.umd.edu's password:"
send "123\r"
expect eof
EOD
exit
fi

make >> /Users/amos/$1.log 2>&1
if [ $? -ne 0 ]
then
cp /Users/amos/$1.log /Users/amos/$1_Failed.log
echo "FAILED: make" >> /Users/amos/$1_Failed.log
/usr/bin/expect <<EOD
spawn scp /Users/amos/$1_Failed.log ssh@sauron.cs.umd.edu:VMlogs
expect "ssh@sauron.cs.umd.edu's password:"
send "123\r"
expect eof
EOD
exit
fi

make check >> /Users/amos/$1.log 2>&1
if [ $? -ne 0 ]
then
cp /Users/amos/$1.log /Users/amos/$1_Failed.log
echo "FAILED: make check" >> /Users/amos/$1_Failed.log
/usr/bin/expect <<EOD
spawn scp /Users/amos/$1_Failed.log ssh@sauron.cs.umd.edu:VMlogs
expect "ssh@sauron.cs.umd.edu's password:"
send "123\r"
expect eof
EOD
exit
fi

echo "AMOS" | sudo -S make install >> /Users/amos/$1.log 2>&1
if [ $? -ne 0 ]
then
cp /Users/amos/$1.log /Users/amos/$1_Failed.log
echo "FAILED: make install" >> /Users/amos/$1_Failed.log
/usr/bin/expect <<EOD
spawn scp /Users/amos/$1_Failed.log ssh@sauron.cs.umd.edu:VMlogs
expect "ssh@sauron.cs.umd.edu's password:"
send "123\r"
expect eof
EOD
exit
fi

echo "AMOS" | sudo -S ln -s /usr/local/AMOS/bin/* /usr/local/bin/

export PATH=$PATH:/usr/local/AMOS/bin
cd test/
./test.sh >> /Users/amos/$1.log 2>&1
if [ $? -ne 0 ]
then
cp /Users/amos/$1.log /Users/amos/$1_Failed.log
echo "FAILED: test.sh" >> /Users/amos/$1_Failed.log
/usr/bin/expect <<EOD
spawn scp /Users/amos/$1_Failed.log ssh@sauron.cs.umd.edu:VMlogs
expect "ssh@sauron.cs.umd.edu's password:"
send "123\r"
expect eof
EOD
exit
fi

./test_amosvalidate.sh >> /Users/amos/$1.log 2>&1
if [ $? -ne 0 ]
then
cp /Users/amos/$1.log /Users/amos/$1_Failed.log
echo "FAILED: test_amosvalidate.sh" >> /Users/amos/$1_Failed.log
/usr/bin/expect <<EOD
spawn scp /Users/amos/$1_Failed.log ssh@sauron.cs.umd.edu:VMlogs
expect "ssh@sauron.cs.umd.edu's password:"
send "123\r"
expect eof
EOD
exit
fi

echo "sending log to walnut..."
now=$(date +"%y%m%d")
echo "SUCCESS:" >> /Users/amos/$1.log
/usr/bin/expect <<EOD
spawn scp /Users/amos/$1.log ssh@sauron.cs.umd.edu:VMlogs
expect "ssh@sauron.cs.umd.edu's password:"
send "123\r"
expect eof
EOD
exit



