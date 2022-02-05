package pal

import (
	"github.com/evolbioinf/fasta"
	"io/ioutil"
	"os"
	"strconv"
	"testing"
)

func TestScoreMatrix(t *testing.T) {
	get := make([]string, 0)
	want := make([]string, 0)
	sm := NewScoreMatrix(1, -3)
	get = append(get, sm.String())
	bf := "./doc/blosum62.txt"
	f, err := os.Open(bf)
	if err != nil {
		t.Errorf("can't open %q", bf)
	}
	defer f.Close()
	sm = ReadScoreMatrix(f)
	get = append(get, sm.String())
	for i, _ := range get {
		f := "data/sm" + strconv.Itoa(i+1) + ".txt"
		b, err := ioutil.ReadFile(f)
		if err != nil {
			t.Errorf("can't open %q", f)
		}
		want = append(want, string(b))
	}
	for i, g := range get {
		w := want[i]
		if g != w {
			t.Errorf("get:\n%s\nwant:\n%s", g, w)
		}
	}
}
func TestAlignment(t *testing.T) {
	get := make([]string, 0)
	want := make([]string, 0)
	bf := "./doc/blosum62.txt"
	f, err := os.Open(bf)
	if err != nil {
		t.Errorf("can't open %q", bf)
	}
	defer f.Close()
	sm := ReadScoreMatrix(f)
	q := fasta.NewSequence("Q", []byte("MKFLALF"))
	s := fasta.NewSequence("S", []byte("MKYLILLF"))
	a := NewGlobalAlignment(q, s, sm, -5, -2)
	a.Align()
	get = append(get, a.String()+"\n")
	s1 := "CATTGGGATGATGCACCTATATTGTGAGGTCTGTTACACTG" +
		"TCGTTCCGCAGATCGAGCAATCCCGTATCCTTTACATAT" +
		"TGCCGTGTGGGGTAAGGTGC"
	s2 := "TTGTGAGGTCTGTTACACTGTCCGCAGATCGAGCAATCC" +
		"CGTATCCTTTACATATTGCCGTGTGGGGTAAGGTGCC" +
		"ATTGGGATGATGCACCTATA"
	q = fasta.NewSequence("Q", []byte(s1))
	s = fasta.NewSequence("S", []byte(s2))
	sm = NewScoreMatrix(1, -3)
	o := NewOverlapAlignment(q, s, sm, -5, -2)
	o.Align()
	get = append(get, o.String()+"\n")
	fn := "data/dmAdhAdhdup.fasta"
	f, err = os.Open(fn)
	if err != nil {
		t.Errorf("can't open %q", fn)
	}
	defer f.Close()
	sc := fasta.NewScanner(f)
	sc.ScanSequence()
	q = sc.Sequence()
	fn = "data/dgAdhAdhdup.fasta"
	f, err = os.Open(fn)
	if err != nil {
		t.Errorf("can't open %q", fn)
	}
	defer f.Close()
	sc = fasta.NewScanner(f)
	sc.ScanSequence()
	s = sc.Sequence()
	l := NewLocalAlignment(q, s, sm, -5, -2)
	l.Align()
	str := l.String()
	l.Align()
	str += "\n" + l.String() + "\n"
	get = append(get, str)
	str = a.PrintMatrix('v')
	get = append(get, str)
	str = a.PrintMatrix('t')
	get = append(get, str)
	for i, _ := range get {
		f := "data/a" + strconv.Itoa(i+1) + ".txt"
		b, err := ioutil.ReadFile(f)
		if err != nil {
			t.Errorf("can't open %q", f)
		}
		want = append(want, string(b))
	}
	for i, g := range get {
		w := want[i]
		if g != w {
			t.Errorf("get:\n%s\nwant:\n%s", g, w)
		}
	}
}
