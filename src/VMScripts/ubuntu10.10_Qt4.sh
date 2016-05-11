#!/bin/sh
cd /
cd home/bryanta/
cd amos/

./bootstrap > /home/bryanta/$1.log 2>&1
if [ $? -ne 0 ]
then
cp /home/bryanta/$1.log /home/bryanta/$1_Failed.log
echo "FAILED: ./bootstrap" >> /home/bryanta/$1_Failed.log
/usr/bin/expect <<EOD
spawn scp /home/bryanta/$1_Failed.log ssh@sauron.cs.umd.edu:VMlogs
expect "ssh@sauron.cs.umd.edu's password:"
send "123\r"
expect eof
EOD
echo "1234561" | sudo -S shutdown -h now
fi

./configure --with-qmake-qt4=/usr/bin/qmake-qt4 --prefix=/usr/local/AMOS >> /home/bryanta/$1.log 2>&1
if [ $? -ne 0 ]
then
cp /home/bryanta/$1.log /home/bryanta/$1_Failed.log
echo "FAILED: ./configure" >> /home/bryanta/$1_Failed.log
/usr/bin/expect <<EOD
spawn scp /home/bryanta/$1_Failed.log ssh@sauron.cs.umd.edu:VMlogs
expect "ssh@sauron.cs.umd.edu's password:"
send "123\r"
expect eof
EOD
echo "1234561" | sudo -S shutdown -h now
fi

make >> /home/bryanta/$1.log 2>&1
if [ $? -ne 0 ]
then
cp /home/bryanta/$1.log /home/bryanta/$1_Failed.log
echo "FAILED: make" >> /home/bryanta/$1_Failed.log
/usr/bin/expect <<EOD
spawn scp /home/bryanta/$1_Failed.log ssh@sauron.cs.umd.edu:VMlogs
expect "ssh@sauron.cs.umd.edu's password:"
send "123\r"
expect eof
EOD
echo "1234561" | sudo -S shutdown -h now
fi

make check >> /home/bryanta/$1.log 2>&1
if [ $? -ne 0 ]
then
cp /home/bryanta/$1.log /home/bryanta/$1_Failed.log
echo "FAILED: make check" >> /home/bryanta/$1_Failed.log
/usr/bin/expect <<EOD
spawn scp /home/bryanta/$1_Failed.log ssh@sauron.cs.umd.edu:VMlogs
expect "ssh@sauron.cs.umd.edu's password:"
send "123\r"
expect eof
EOD
echo "1234561" | sudo -S shutdown -h now
fi

echo "1234561" | sudo -S make install >> /home/bryanta/$1.log 2>&1
if [ $? -ne 0 ]
then
cp /home/bryanta/$1.log /home/bryanta/$1_Failed.log
echo "FAILED: make install" >> /home/bryanta/$1_Failed.log
/usr/bin/expect <<EOD
spawn scp /home/bryanta/$1_Failed.log ssh@sauron.cs.umd.edu:VMlogs
expect "ssh@sauron.cs.umd.edu's password:"
send "123\r"
expect eof
EOD
echo "1234561" | sudo -S shutdown -h now
fi

echo "1234561" | sudo -S ln -s /usr/local/AMOS/bin/* /usr/local/bin/
export PATH=$PATH:/usr/local/AMOS/bin
cd test/
./test.sh >> /home/bryanta/$1.log 2>&1
if [ $? -ne 0 ]
then
cp /home/bryanta/$1.log /home/bryanta/$1_Failed.log
echo "FAILED: tesh.sh" >> /home/bryanta/$1_Failed.log
/usr/bin/expect <<EOD
spawn scp /home/bryanta/$1_Failed.log ssh@sauron.cs.umd.edu:VMlogs
expect "ssh@sauron.cs.umd.edu's password:"
send "123\r"
expect eof
EOD
echo "1234561" | sudo -S shutdown -h now
sleep 180
fi

./test_amosvalidate.sh >> /home/bryanta/$1.log 2>&1
if [ $? -ne 0 ]
then
cp /home/bryanta/$1.log /home/bryanta/$1_Failed.log
echo "FAILED: test_amosvalidate.sh" >> /home/bryanta/$1_Failed.log
/usr/bin/expect <<EOD
spawn scp /home/bryanta/$1_Failed.log ssh@sauron.cs.umd.edu:VMlogs
expect "ssh@sauron.cs.umd.edu's password:"
send "123\r"
expect eof
EOD
echo "1234561" | sudo -S shutdown -h now
sleep 180
fi

echo "sending log to walnut..."
now=$(date +"%y%m%d")
echo "SUCCESS:" >> /home/bryanta/$1.log
/usr/bin/expect <<EOD
spawn scp /home/bryanta/$1.log ssh@sauron.cs.umd.edu:VMlogs
expect "ssh@sauron.cs.umd.edu's password:"
send "123\r"
expect eof
EOD

echo "deleting log..."
rm /home/bryanta/$1.log

echo "shutting down..."
echo "1234561" | sudo -S shutdown -h now
