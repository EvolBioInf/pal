//  Package pal implements data structures, methods, and functions  for computing and printing optimal pairwise alignments.
package pal

import (
	"bufio"
	"bytes"
	"fmt"
	"github.com/evolbioinf/fasta"
	"io"
	"log"
	"math"
	"sort"
	"strconv"
	"text/tabwriter"
)

//  A ScoreMatrix stores the scores of residue pairs.
type ScoreMatrix struct {
	res string
	mat [][]float64
	dic map[byte]int
}
type alignment struct {
	q, s             *fasta.Sequence
	m                *ScoreMatrix
	gapO, gapE       float64
	ql, sl           int
	p                [][]Cell
	qa, sa           []byte
	score            float64
	gaps, mismatches int
	ll               int
	qs               int
	ss               int
	count            int
	coords           []coordinate
}
type Cell struct {
	E, F, G, V float64
	visited    bool
}

// GlobalAlignment holds the global alignment of two sequences.
type GlobalAlignment struct {
	alignment
}

// OverlapAlignment holds the overlap alignment of two sequences.
type OverlapAlignment struct {
	alignment
}

// LocalAlignment holds the local alignment of two sequences.
type LocalAlignment struct {
	alignment
}
type coordinate struct {
	i, j int
	s    float64
}
type coordSlice []coordinate

func (s *ScoreMatrix) setScore(r1, r2 byte, sc float64) {
	if r1 >= 97 {
		r1 -= 32
	}
	if r2 >= 97 {
		r2 -= 32
	}
	i1, ok1 := s.dic[r1]
	i2, ok2 := s.dic[r2]
	if !ok1 || !ok2 {
		log.Fatalf("can't set score for "+
			"residues (%c, %c)", r1, r2)
	}
	s.mat[i1][i2] = sc
}

// The method Score takes two characters as arguments and returns their score. If one of the characters is not a printing ASCII character, it exits with message.
func (s *ScoreMatrix) Score(r1, r2 byte) float64 {
	if r1 >= 97 {
		r1 -= 32
	}
	if r2 >= 97 {
		r2 -= 32
	}
	i1, ok1 := s.dic[r1]
	i2, ok2 := s.dic[r2]
	if !ok1 || !ok2 {
		log.Fatalf("can't score (%c, %c)", r1, r2)
	}
	return s.mat[i1][i2]
}

// String converts the ScoreMatrix to a printable string.
func (s *ScoreMatrix) String() string {
	b := make([]byte, 0)
	buf := bytes.NewBuffer(b)
	w := tabwriter.NewWriter(buf, 1, 1, 1, ' ', 0)
	for _, r := range s.res {
		fmt.Fprintf(w, "\t%2c", r)
	}
	fmt.Fprintf(w, "\n")
	n := len(s.res)
	for i := 0; i < n; i++ {
		fmt.Fprintf(w, "%c", s.res[i])
		for j := 0; j < n; j++ {
			fmt.Fprintf(w, "\t%2.3g", s.mat[i][j])
		}
		if i < n-1 {
			fmt.Fprintf(w, "\n")
		}
	}
	w.Flush()
	return buf.String()
}
func (a *alignment) new(q, s *fasta.Sequence,
	sm *ScoreMatrix,
	gapO, gapE float64) {
	a.q = q
	a.s = s
	a.m = sm
	a.gapO = gapO
	a.gapE = gapE
	a.ql = len(q.Data())
	a.sl = len(s.Data())
	m := len(a.q.Data()) + 1
	n := len(a.s.Data()) + 1
	a.p = make([][]Cell, m)
	for i := 0; i < m; i++ {
		a.p[i] = make([]Cell, n)
	}

	a.qa = make([]byte, 0)
	a.sa = make([]byte, 0)
	a.ll = fasta.DefaultLineLength
	a.coords = make([]coordinate, 0)
}

//  The Method ProgrammingMatrix returns the (m+1) x (n+1) matrix  used for calculating an alignment of two sequences of lengths m and  n.
func (a *alignment) ProgrammingMatrix() [][]Cell {
	return a.p
}

// RawAlignment returns the aligned query and subject sequences.
func (a *alignment) RawAlignment() (q, s []byte) {
	return a.qa, a.sa
}

// Method String returns the alignment as a string ready for pretty printing.
func (a *alignment) String() string {
	var buf []byte
	buffer := bytes.NewBuffer(buf)
	w := new(tabwriter.Writer)
	w.Init(buffer, 1, 0, 1, ' ', 0)
	h := a.q.Header()
	fmt.Fprintf(w, "Query\t%s\t(%d residue", h, a.ql)
	if a.ql != 1 {
		fmt.Fprint(w, "s")
	}
	fmt.Fprint(w, ")\n")
	h = a.s.Header()
	fmt.Fprintf(w, "Subject\t%s\t(%d residue", h, a.sl)
	if a.sl != 1 {
		fmt.Fprint(w, "s")
	}
	fmt.Fprint(w, ")\n")
	fmt.Fprintf(w, "Score\t%g\n", a.score)
	er := a.gaps + a.mismatches
	fmt.Fprintf(w, "Error")
	if er != 1 {
		fmt.Fprintf(w, "s")
	}
	fmt.Fprintf(w, "\t%d", er)
	fmt.Fprintf(w, " (%d gap", a.gaps)
	if a.gaps != 1 {
		fmt.Fprint(w, "s")
	}
	fmt.Fprintf(w, ", %d mismatch", a.mismatches)
	if a.mismatches != 1 {
		fmt.Fprint(w, "es")
	}
	fmt.Fprint(w, ")\n")

	w.Flush()
	end := 0
	qs := a.qs
	matches := make([]rune, a.ll)
	ss := a.ss
	for i := 0; i < len(a.qa); i += a.ll {
		if i+a.ll < len(a.qa) {
			end = i + a.ll
		} else {
			end = len(a.qa)
		}
		data := a.qa[i:end]
		nr := len(data) - bytes.Count(data, []byte("-"))
		l := qs
		if nr > 0 {
			l++
		}
		fmt.Fprintf(w, "\n\nQuery\t%d\t%s\t%d\n", l, data, qs+nr)
		qs += nr
		k := 0
		for j := i; j < end; j++ {
			m := ' '
			qc := a.qa[j]
			sc := a.sa[j]
			if qc != '-' && sc != '-' {
				if qc == sc {
					m = '|'
				} else if a.m.Score(qc, sc) > 0 {
					m = ':'
				}
			}
			matches[k] = m
			k++
		}
		fmt.Fprintf(w, "\t\t%s\n", string(matches[:k]))
		data = a.sa[i:end]
		nr = len(data) - bytes.Count(data, []byte("-"))
		l = ss
		if nr > 0 {
			l++
		}
		fmt.Fprintf(w, "Subject\t%d\t%s\t%d\n", l, data, ss+nr)
		ss += nr
	}
	w.Flush()
	buffer.Write([]byte("//"))
	return buffer.String()
}
func (a *alignment) reverse() {
	revStr(a.qa)
	revStr(a.sa)
}
func (a *alignment) errors() {
	a.mismatches = 0
	a.gaps = 0
	if len(a.qa) != len(a.sa) {
		log.Fatal("aligned sequences don't " +
			"have same length")
	}
	for i, r1 := range a.qa {
		r2 := a.sa[i]
		if r1 == '-' || r2 == '-' {
			a.gaps++
		} else if r1 != r2 {
			a.mismatches++
		}
	}
}

// SetLineLength sets the length of alignment lines in String.
func (a *alignment) SetLineLength(ll int) {
	if ll > 0 {
		a.ll = ll
	}
}

// Score returns the score of the alignment.
func (a *alignment) Score() float64 {
	return a.score
}

// Mismatches returns the number of mismatches of the alignment.
func (a *alignment) Mismatches() int {
	return a.mismatches
}

// Gaps returns the number of gap characters in the alignment.
func (a *alignment) Gaps() int {
	return a.gaps
}

// SetSubjectStart sets the start of the subject sequence, which is useful in global/local alignment.
func (a *alignment) SetSubjectStart(s int) {
	a.ss = s
}

// SetSubjectLength sets the length of the subject sequence, which is useful in global/local alignment.
func (a *alignment) SetSubjectLength(l int) {
	a.sl = l
}

// Method Align computes the alignment.
func (a *GlobalAlignment) Align() {
	q := a.q.Data()
	s := a.s.Data()
	m := len(q)
	n := len(s)
	p := a.p
	for i := 1; i <= m; i++ {
		p[i][0].E = a.gapO + float64(i-1)*a.gapE
		p[i][0].V = p[i][0].E
		p[i][0].F = -math.MaxFloat64
		p[i][0].G = -math.MaxFloat64
	}
	for j := 1; j <= n; j++ {
		p[0][j].F = a.gapO + float64(j-1)*a.gapE
		p[0][j].V = p[0][j].F
		p[0][j].E = -math.MaxFloat64
		p[0][j].G = -math.MaxFloat64
	}
	for i := 1; i <= m; i++ {
		for j := 1; j <= n; j++ {
			p[i][j].E = math.Max(p[i-1][j].E+a.gapE, p[i-1][j].V+a.gapO)
			p[i][j].F = math.Max(p[i][j-1].F+a.gapE, p[i][j-1].V+a.gapO)
			p[i][j].G = p[i-1][j-1].V + a.m.Score(q[i-1], s[j-1])
			p[i][j].V = math.Max(math.Max(p[i][j].E, p[i][j].F), p[i][j].G)
		}
	}
	a.score = a.p[m][n].V
	i := m
	j := n
	for i > 0 || j > 0 {
		if p[i][j].V == p[i][j].E {
			a.qa = append(a.qa, q[i-1])
			a.sa = append(a.sa, '-')
			i--
		} else if p[i][j].V == p[i][j].F {
			a.qa = append(a.qa, '-')
			a.sa = append(a.sa, s[j-1])
			j--
		} else {
			a.qa = append(a.qa, q[i-1])
			a.sa = append(a.sa, s[j-1])
			i--
			j--
		}
	}
	a.reverse()
	a.errors()
}
func (a *OverlapAlignment) Align() {
	q := a.q.Data()
	s := a.s.Data()
	m := len(q)
	n := len(s)
	p := a.p
	for i := 1; i <= m; i++ {
		for j := 1; j <= n; j++ {
			p[i][j].E = math.Max(p[i-1][j].E+a.gapE, p[i-1][j].V+a.gapO)
			p[i][j].F = math.Max(p[i][j-1].F+a.gapE, p[i][j-1].V+a.gapO)
			p[i][j].G = p[i-1][j-1].V + a.m.Score(q[i-1], s[j-1])
			p[i][j].V = math.Max(math.Max(p[i][j].E, p[i][j].F), p[i][j].G)
		}
	}
	a.score = a.p[m][n].V
	max := math.Inf(-1)
	j := 0
	i := m
	for k := 0; k <= n; k++ {
		if max < p[i][k].V {
			max = p[i][k].V
			j = k
		}
	}
	for k := 0; k <= m; k++ {
		if max < p[k][n].V {
			max = p[k][n].V
			i = k
			j = n
		}
	}
	a.score = p[i][j].V
	for k := m; k > i; k-- {
		a.qa = append(a.qa, q[k-1])
		a.sa = append(a.sa, '-')
	}
	for k := n; k > j; k-- {
		a.qa = append(a.qa, '-')
		a.sa = append(a.sa, s[k-1])
	}
	for i > 0 && j > 0 {
		if p[i][j].V == p[i][j].E {
			a.qa = append(a.qa, q[i-1])
			a.sa = append(a.sa, '-')
			i--
		} else if p[i][j].V == p[i][j].F {
			a.qa = append(a.qa, '-')
			a.sa = append(a.sa, s[j-1])
			j--
		} else {
			a.qa = append(a.qa, q[i-1])
			a.sa = append(a.sa, s[j-1])
			i--
			j--
		}
	}
	for k := i; k > 0; k-- {
		a.sa = append(a.sa, '-')
		a.qa = append(a.qa, q[k-1])
	}
	for k := j; k > 0; k-- {
		a.sa = append(a.sa, s[k-1])
		a.qa = append(a.qa, '-')
	}
	a.reverse()
	a.errors()
}

// TrimQuery removes gaps flanking the query.
func (a *OverlapAlignment) TrimQuery() {
	if len(a.qa) != len(a.sa) {
		log.Fatal("can't trim alignments " +
			"of unequal length")
	}
	for a.qa[0] == '-' {
		a.qa = a.qa[1:]
		a.sa = a.sa[1:]
		a.ss++
		a.gaps--
	}
	for a.qa[len(a.qa)-1] == '-' {
		a.qa = a.qa[:len(a.qa)-1]
		a.sa = a.sa[:len(a.sa)-1]
		a.gaps--
	}
}

// Align computes the next alignment in the dynamic programming matrix. It returns false if none was found.
func (a *LocalAlignment) Align() bool {
	q := a.q.Data()
	s := a.s.Data()
	m := len(q)
	n := len(s)
	p := a.p
	a.count++
	if a.count == 1 {
		for i := 1; i <= m; i++ {
			for j := 1; j <= n; j++ {
				p[i][j].E = math.Max(p[i-1][j].E+a.gapE, p[i-1][j].V+a.gapO)
				p[i][j].F = math.Max(p[i][j-1].F+a.gapE, p[i][j-1].V+a.gapO)
				p[i][j].G = p[i-1][j-1].V + a.m.Score(q[i-1], s[j-1])
				p[i][j].V = math.Max(math.Max(p[i][j].E, p[i][j].F), p[i][j].G)
				if p[i][j].V < 0 {
					p[i][j].V = 0
				}
			}
		}
		var c coordinate
		c.s = -math.MaxFloat64
		for i := 1; i <= m; i++ {
			for j := 1; j <= n; j++ {
				if c.s < p[i][j].V {
					c.s = p[i][j].V
					c.i = i
					c.j = j
				}
			}
		}
		a.coords = append(a.coords, c)
	} else if a.count == 2 {
		var c coordinate
		for i := 1; i <= m; i++ {
			for j := 1; j <= n; j++ {
				if !p[i][j].visited {
					c.i = i
					c.j = j
					c.s = p[i][j].V
					a.coords = append(a.coords, c)
				}
			}
		}
		sort.Sort(coordSlice(a.coords))
	}
	found := false
	for k, c := range a.coords {
		i := c.i
		j := c.j
		found = true
		a.score = p[i][j].V
		for p[i][j].V > 0 {
			p[i][j].visited = true
			if p[i][j].V == p[i][j].E {
				a.qa = append(a.qa, q[i-1])
				a.sa = append(a.sa, '-')
				i--
			} else if p[i][j].V == p[i][j].F {
				a.qa = append(a.qa, '-')
				a.sa = append(a.sa, s[j-1])
				j--
			} else {
				a.qa = append(a.qa, q[i-1])
				a.sa = append(a.sa, s[j-1])
				i--
				j--
			}
			if p[i][j].visited {
				found = false
				a.qa = a.qa[:0]
				a.sa = a.sa[:0]
				break
			}
		}
		if found {
			a.qs = i
			a.ss = j
			a.reverse()
			a.errors()

		}
		if found {
			a.coords = a.coords[k+1:]
			break
		}
	}
	return found
}
func (c coordSlice) Len() int {
	return len(c)
}
func (c coordSlice) Swap(i, j int) {
	c[i], c[j] = c[j], c[i]
}
func (c coordSlice) Less(i, j int) bool {
	return c[i].s > c[j].s
}

// Function NewScoreMatrix generates a new score matrix from a match and a mismatch score.
func NewScoreMatrix(match, mismatch float64) *ScoreMatrix {
	sm := new(ScoreMatrix)
	sm.dic = make(map[byte]int)
	sm.res = "ARNDCQEGHIJLKMFPSTWYVBZX*U"
	for i, r := range sm.res {
		sm.dic[byte(r)] = i
	}
	n := len(sm.res)
	sm.mat = make([][]float64, n)
	for i := 0; i < n; i++ {
		sm.mat[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		sm.mat[i][i] = match
	}
	for i := 0; i < n-1; i++ {
		for j := i + 1; j < n; j++ {
			sm.mat[i][j] = mismatch
			sm.mat[j][i] = mismatch
		}
	}
	return sm
}

// Function ReadScores reads the entries of a ScoreMatrix from an io.Reader.
func ReadScoreMatrix(r io.Reader) *ScoreMatrix {
	s := NewScoreMatrix(1, -1)
	first := true
	var res [][]byte
	sc := bufio.NewScanner(r)
	for sc.Scan() {
		b := sc.Bytes()
		if b[0] != '#' {
			if first {
				res = bytes.Fields(b)
				first = false
			} else {
				entries := bytes.Fields(b)
				for i := 1; i < len(entries); i++ {
					r1 := entries[0][0]
					r2 := res[i-1][0]
					score, err := strconv.ParseFloat(string(entries[i]), 64)
					if err != nil {
						log.Fatalf("couldn't parse %q\n", entries[i])
					}
					s.setScore(r1, r2, score)
				}
			}
		}
	}
	return s
}
func revStr(s []byte) {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
}

// NewGlobalAlignment constructs a new global alignment from a query and a subject sequence, and a substitution matrix. To compute the actual alignment, call the method Align.
func NewGlobalAlignment(q, s *fasta.Sequence,
	sm *ScoreMatrix,
	gapO, gapE float64) *GlobalAlignment {
	ga := new(GlobalAlignment)
	ga.new(q, s, sm, gapO, gapE)
	return ga
}

// NewOverlapAlignment constructs an overlap alignment from a query and a subject sequence, and a substitution matrix. To compute the actual alignment, call the method Align.
func NewOverlapAlignment(q, s *fasta.Sequence,
	sm *ScoreMatrix,
	gapO, gapE float64) *OverlapAlignment {
	oa := new(OverlapAlignment)
	oa.new(q, s, sm, gapO, gapE)
	return oa
}

// NewLocalAlignment constructs a local alignment from a query and a subject sequence, and a substitution matrix. To compute the actual alignment, call the method Align.
func NewLocalAlignment(q, s *fasta.Sequence,
	sm *ScoreMatrix,
	gapO, gapE float64) *LocalAlignment {
	la := new(LocalAlignment)
	la.new(q, s, sm, gapO, gapE)
	return la
}
