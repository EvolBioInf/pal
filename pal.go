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
	q, s       *fasta.Sequence
	m          *ScoreMatrix
	gapO, gapE float64
	p          [][]cell
	qa, sa     []byte
	score      float64
	ll         int
	qs         int
	ss         int
	count      int
	coords     []coordinate
}
type cell struct {
	e, f, g, v float64
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
	m := len(a.q.Data()) + 1
	n := len(a.s.Data()) + 1
	a.p = make([][]cell, m)
	for i := 0; i < m; i++ {
		a.p[i] = make([]cell, n)
	}

	a.qa = make([]byte, 0)
	a.sa = make([]byte, 0)
	a.ll = fasta.DefaultLineLength
	a.coords = make([]coordinate, 0)
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
	l := len(a.q.Data())
	fmt.Fprintf(w, "Query\t%s\t(%d residues)\n", h, l)
	h = a.s.Header()
	l = len(a.s.Data())
	fmt.Fprintf(w, "Subject\t%s\t(%d residues)\n", h, l)
	fmt.Fprintf(w, "Score\t%g\n", a.score)
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

// Method Align computes the alignment.  We declare a set of
func (a *GlobalAlignment) Align() {
	q := a.q.Data()
	s := a.s.Data()
	m := len(q)
	n := len(s)
	p := a.p
	for i := 1; i <= m; i++ {
		p[i][0].e = a.gapO + float64(i)*a.gapE
		p[i][0].v = p[i][0].e
		p[i][0].f = -math.MaxFloat64
		p[i][0].g = -math.MaxFloat64
	}
	for j := 1; j <= n; j++ {
		p[0][j].f = a.gapO + float64(j)*a.gapE
		p[0][j].v = p[0][j].f
		p[0][j].e = -math.MaxFloat64
		p[0][j].g = -math.MaxFloat64
	}
	for i := 1; i <= m; i++ {
		for j := 1; j <= n; j++ {
			p[i][j].e = math.Max(p[i-1][j].e, p[i-1][j].v+a.gapO) + a.gapE
			p[i][j].f = math.Max(p[i][j-1].f, p[i][j-1].v+a.gapO) + a.gapE
			p[i][j].g = p[i-1][j-1].v + a.m.Score(q[i-1], s[j-1])
			p[i][j].v = math.Max(math.Max(p[i][j].e, p[i][j].f), p[i][j].g)
		}
	}
	a.score = a.p[m][n].v
	i := m
	j := n
	for i > 0 || j > 0 {
		if p[i][j].v == p[i][j].e {
			a.qa = append(a.qa, q[i-1])
			a.sa = append(a.sa, '-')
			i--
		} else if p[i][j].v == p[i][j].f {
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
	reverse(a.qa)
	reverse(a.sa)
}
func (a *OverlapAlignment) Align() {
	q := a.q.Data()
	s := a.s.Data()
	m := len(q)
	n := len(s)
	p := a.p
	for i := 1; i <= m; i++ {
		for j := 1; j <= n; j++ {
			p[i][j].e = math.Max(p[i-1][j].e, p[i-1][j].v+a.gapO) + a.gapE
			p[i][j].f = math.Max(p[i][j-1].f, p[i][j-1].v+a.gapO) + a.gapE
			p[i][j].g = p[i-1][j-1].v + a.m.Score(q[i-1], s[j-1])
			p[i][j].v = math.Max(math.Max(p[i][j].e, p[i][j].f), p[i][j].g)
		}
	}
	a.score = a.p[m][n].v
	max := math.Inf(-1)
	j := 0
	i := m
	for k := 0; k <= n; k++ {
		if max < p[i][k].v {
			max = p[i][k].v
			j = k
		}
	}
	for k := 0; k <= m; k++ {
		if max < p[k][n].v {
			max = p[k][n].v
			i = k
			j = n
		}
	}
	a.score = p[i][j].v
	for k := m; k > i; k-- {
		a.qa = append(a.qa, q[k-1])
		a.sa = append(a.sa, '-')
	}
	for k := n; k > j; k-- {
		a.qa = append(a.qa, '-')
		a.sa = append(a.sa, s[k-1])
	}
	for i > 0 && j > 0 {
		if p[i][j].v == p[i][j].e {
			a.qa = append(a.qa, q[i-1])
			a.sa = append(a.sa, '-')
			i--
		} else if p[i][j].v == p[i][j].f {
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
				p[i][j].e = math.Max(p[i-1][j].e, p[i-1][j].v+a.gapO) + a.gapE
				p[i][j].f = math.Max(p[i][j-1].f, p[i][j-1].v+a.gapO) + a.gapE
				p[i][j].g = p[i-1][j-1].v + a.m.Score(q[i-1], s[j-1])
				p[i][j].v = math.Max(math.Max(p[i][j].e, p[i][j].f), p[i][j].g)
				if p[i][j].v < 0 {
					p[i][j].v = 0
				}
			}
		}
		var c coordinate
		c.s = -math.MaxFloat64
		for i := 1; i <= m; i++ {
			for j := 1; j <= n; j++ {
				if c.s < p[i][j].v {
					c.s = p[i][j].v
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
					c.s = p[i][j].v
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
		a.score = p[i][j].v
		for p[i][j].v > 0 {
			p[i][j].visited = true
			if p[i][j].v == p[i][j].e {
				a.qa = append(a.qa, q[i-1])
				a.sa = append(a.sa, '-')
				i--
			} else if p[i][j].v == p[i][j].f {
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
	sm.res = "ARNDCQEGHILKMFPSTWYVBZX*U"
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

// NewGlobalAlignment constructs a new global alignment from a query and a subject sequence, and a substitution matrix. To compute the actual alignment, invoke the method Align.
func NewGlobalAlignment(q, s *fasta.Sequence,
	sm *ScoreMatrix,
	gapO, gapE float64) *GlobalAlignment {
	ga := new(GlobalAlignment)
	ga.new(q, s, sm, gapO, gapE)
	return ga
}
func reverse(s []byte) {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
}

// NewOverlapAlignment constructs an overlap alignment from a query and a subject sequence, and a substitution matrix. To compute the actual alignment, invoke the method Align.
func NewOverlapAlignment(q, s *fasta.Sequence,
	sm *ScoreMatrix,
	gapO, gapE float64) *OverlapAlignment {
	oa := new(OverlapAlignment)
	oa.new(q, s, sm, gapO, gapE)
	return oa
}

// NewLocalAlignment constructs a local alignment from a query and a subject sequence, and a substitution matrix. To compute the actual alignment, invoke the method Align.
func NewLocalAlignment(q, s *fasta.Sequence,
	sm *ScoreMatrix,
	gapO, gapE float64) *LocalAlignment {
	la := new(LocalAlignment)
	la.new(q, s, sm, gapO, gapE)
	return la
}
