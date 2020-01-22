#!bin/bash
cd /home/opsrv

# download the latest release
LATEST_RELEASE="$(curl -s https://github.com/saezlab/pypath/releases/latest \
| grep "tag" \
| cut -d \" -f 2 \
| cut -d "/" -f 8)"

curl -sOL "https://github.com/saezlab/pypath/archive/"$LATEST_RELEASE'.tar.gz'

pip install $LATEST_RELEASE.tar.gz
tar -xvf $LATEST_RELEASE.tar.gz
cd pypath-$LATEST_RELEASE/

ps -ef | grep "python ./server-runner-new.py" | grep -v grep | awk '{print $2}' | xargs kill
nohup python server-runner-new.py &
