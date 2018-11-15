
TESTS := test_config.py test_discrete_legendre.py

.PHONY: test

test:
	for t in $(TESTS); do echo "== $$t ==" && (python3 $$t || exit 1); done

