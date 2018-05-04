
# apidoc/ generate by
# sphinx-apidoc -F -o apidoc ../lss-ps/py/lssps

.PHONY: apidoc open

apidoc:
	cd apidoc && make html

open:
	open apidoc/_build/html/index.html
