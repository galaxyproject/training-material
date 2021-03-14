# if a new ubuntu image add some basics first
# apt update ; apt install -y python3-dev python3-venv python3-wheel nano curl wget git python3-setuptools
# adjust to suit your needs. GALDIR could be an existing dev directory, and the curl line could be commented out to save time
# https://github.com/galaxyproject/galaxy/archive/release_21.01.zip
GALDIR="galaxy-central"
GALURL="https://github.com/galaxyproject/galaxy/archive/release_21.01.tar.gz"
# https://github.com/galaxyproject/galaxy/archive/dev.tar.gz if you want dev
PDIR="planemo"
CDIR=`pwd`
git clone --recursive https://github.com/fubar2/planemo.git $PDIR
rm -rf $PDIR/docs
mkdir -p $GALDIR
curl -L -s $GALURL | tar xzf - --strip-components=1 -C $GALDIR
cp $PDIR/planemo_ext/welcome.html $GALDIR/static/welcome.html
cp $PDIR/planemo_ext/welcome.html $GALDIR/static/welcome.html.sample
sed -i '/-l|-list|--list)/i\\t --dev-wheels|-dev-wheels)\n\t\t shift\n\t\t ;;\n' $GALDIR/run_tests.sh
# planemo will not run as a tool successfully without this - something fishy with recent changes to planemo and run_tests.sh. I dunno
# https://github.com/galaxyproject/planemo/issues/1148
sed 's/#sanitize_all_html\: true/sanitize_all_html\: false/g' $GALDIR/config/galaxy.yml.sample > $GALDIR/config/galaxy.yml
# need this to see html in collections like planemo tests
mkdir -p $PDIR/mytools
cd $PDIR
python3 -m venv .venv
. .venv/bin/activate
python3 setup.py build
python3 setup.py install
cd $CDIR
planemo conda_init --conda_prefix $PDIR/con
planemo tool_factory --galaxy_root $GALDIR --port 9090 --host 0.0.0.0 --conda_prefix $PDIR/con
# use planemo tool_factory --galaxy_root galaxy-central --port 9090 --host 0.0.0.0 --conda_prefix planemo/con
# after activating the venv as above to restart planemo next time without all the downloading
# ALL YOUR WORK WILL BE GONE unless you explicitly exported your ToolFactory jobs as histories or as workflows.


#git clone --recursive https://github.com/fubar2/planemo.git
#pip install ./planemo
#planemo tool_factory --port 9090 --host 0.0.0.0 --install_galaxy --conda_auto_install

#google-api-core 1.25.0 requires google-auth<2.0dev,>=1.21.1, but you have google-auth 1.18.0 which is incompatible
# using above in a fresh venv
# (.venv) 1 ross@rosspn50:/tmp/foo$ diff  ../zot/galaxy-central/run_tests.sh galaxy-central/run_tests.sh
# 354a355,357
# >       --dev-wheels|-dev-wheels)
# >           shift
# >           ;;
