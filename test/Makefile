
TESTS := test_config.py 

.PHONY: test

test:
	for t in $(TESTS); do (python3 $$t || exit 1); done

