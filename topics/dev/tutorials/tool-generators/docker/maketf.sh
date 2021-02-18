# if a new ubuntu image, will need a port mapped and add some basics first
# apt update ; apt install -y python3-dev python3-venv python3-wheel nano curl wget git python3-setuptools
TARGDIR="galaxy-central"
PDIR="planemo"
git clone --recursive https://github.com/fubar2/planemo.git $PDIR
rm -rf $PDIR/docs
mkdir -p $TARGDIR
curl -L -s https://github.com/galaxyproject/galaxy/archive/dev.tar.gz | tar xzf - --strip-components=1 -C $TARGDIR
cp $PDIR/planemo_ext/welcome.html $TARGDIR/static/
mkdir -p $PDIR/mytools
cd $PDIR
python3 -m venv .venv
. .venv/bin/activate
python3 setup.py build
python3 setup.py install
planemo conda_init --conda_prefix ./con
planemo tool_factory --galaxy_root $TARGDIR --port 9090 --host 0.0.0.0 --conda_prefix $PDIR/con

# planemo tool_factory --galaxy_root ./galaxy-central --port 8080 --host 0.0.0.0 --conda_prefix /planemo/con
# host is needed to get -p 9090:9090 to work in docker. Default 127.0.0.1 doesn't redirect :(ls -l /tmp
# echo "Starting first run. This takes ages and includes building the Galaxy client. Be patient. Do something else for 20 minutes"
# planemo tool_factory --galaxy_root ./galaxy-central --port 9090 --host 0.0.0.0 --conda_prefix ./planemo/con --conda_auto_init --conda_auto_install
