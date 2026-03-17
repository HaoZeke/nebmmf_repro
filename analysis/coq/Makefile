COQC = coqc
COQFLAGS = -R . OCI_NEB

VFILES = ThresholdBounds.v PenaltyConsistency.v TriggerRedundancy.v MonotonicProgress.v AlgorithmInvariants.v HouseholderEquivalence.v SaddleStability.v
VOFILES = $(VFILES:.v=.vo)

all: $(VOFILES)

%.vo: %.v
	$(COQC) $(COQFLAGS) $<

clean:
	rm -f *.vo *.vok *.vos *.glob .*.aux

.PHONY: all clean
