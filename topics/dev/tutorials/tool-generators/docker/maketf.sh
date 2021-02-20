# if a new ubuntu image add some basics first
# apt update ; apt install -y python3-dev python3-venv python3-wheel nano curl wget git python3-setuptools
# adjust to suit your needs. GALDIR could be an existing dev directory, and the curl line could be commented out to save time
GALDIR="galaxy-central"
PDIR="planemo"
CDIR=`pwd`
git clone --recursive https://github.com/fubar2/planemo.git $PDIR
rm -rf $PDIR/docs
mkdir -p $GALDIR
curl -L -s https://github.com/galaxyproject/galaxy/archive/dev.tar.gz | tar xzf - --strip-components=1 -C $GALDIR
cp $PDIR/planemo_ext/welcome.html $GALDIR/static/welcome.html
cp $PDIR/planemo_ext/welcome.html $GALDIR/static/welcome.html.sample
mkdir -p $PDIR/mytools
cd $PDIR
python3 -m venv .venv
. .venv/bin/activate
python3 setup.py build
python3 setup.py install
cd $CDIR
planemo conda_init --conda_prefix $PDIR/con
planemo tool_factory --galaxy_root $GALDIR --port 8081 --host 0.0.0.0 --conda_prefix $PDIR/con
