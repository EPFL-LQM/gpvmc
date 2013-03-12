ALL: utils

utils: vmc.py stagflux.py setup.py
	python setup.py sdist

install: utils
	python setup.py install --user
