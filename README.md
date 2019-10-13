# Nucleolus

version 1.2 16/08/2019

Various algorithms finding and verifying the nucleolus of cooperative games.

This is a prototype code that has not been developed for distribution.

## Finding the nucleolus

### Algorithms

For each algorithm covered there is an executable and the source as well.

 * BNF: Benedek et al. (2018) - Finding and verifying the nucleolus of cooperative games (forthcoming)
 * DK: [Derks and Kuipers (1997) - Implementing the simplex method for computing the prenucleolus of transferable utility games](https://www.researchgate.net/publication/265080719)
 * SP: the primal sequence of [Solymosi (1993) - On computing the nucleolus of cooperative games](https://www.researchgate.net/publication/318017147)
 * SD: the dual sequence of [Solymosi (1993) - On computing the nucleolus of cooperative games](https://www.researchgate.net/publication/318017147)
 * PD: primal-dual algorithm from Benedek et al. (2019) - Trade-offs in the computation of the nucleolus (forthcoming)
 * DP: dual-primal algorithm from Benedek et al. (2019) - Trade-offs in the computation of the nucleolus (forthcoming)
 * PRA: [Potters et al. (1996) - Computing the Nucleolus by Solving a Prolonged Simplex Algorithm] (https://pubsonline.informs.org/doi/abs/10.1287/moor.21.3.757)
 

### Executables

All executables are obtained by compiling the corresponding source (see below) using Microsoft Visual Studio Enterprise 2015 (version 14.0.25431.01 Update 3) on OS Microsoft Windows 10 Enterprise (version 10.0.17134 Build 17134). NOTE THAT ALL EXECUTABLES ARE USING CPLEX (version 12.7.0)

 * BNF.exe
 * DK.exe
 * SP.exe
 * SD.exe
 * PD.exe
 * DP.exe
 * PRA.exe

### Sources

Each of the codes below requires common.cpp, common.h, gen_game.cpp and gen_game.h
 * BNF.cpp, BNF.h
 * DK.cpp, DK.h
 * SP.cpp, SP.h
 * SD.cpp, SD.h
 * PD.cpp, PD.h
 * DP.cpp, DP.h
 * PRA.cpp, PRA.h

### Inputs

input.txt: containing variables
```
n, type, seed, disp, memo, nlsu
```
separated by linebreaks
 * n: integer, the number of players
 * type: integer between 0 and 7; type=0 reads the game from "v.txt" (see below), type>0 generates certain types of games (see below)
 * seed: integer, seed for randomization when generating certain types of games
 * disp: boolean, 1 to display information while running (default 0)
 * memo: boolean, 1 to switch to memory-saving implementation (default 0); as a rough guide, using 16 Gb memory only memo=1 works with ~n>28 for most games (not applicable for PRA)
 * nlsu: boolean, 1 to switch off linear speed-up, that removes all coalitions in the linear span of the settled coalitions (not applicable for PRA)

v.txt: being read only if type=0
 * containing coalitional values of a game separated by linebreaks
 * for an n-player game it should contain 2^n-1 lines
 * the coalitional value in the i-th row (for i between 1 and 2^n-1) corresponds to the coalition with n-dim. characteristic vector as the binary number converted from the decimal number i
 * for example if n=5, then the file should contain 31 lines
 * the first line should contain the value of player 1
 * the 7-th line should contain the value of coalition containing players 1, 2 and 3
 * the 19-th line should contain the value of coalition containing players 1, 2 and 5

types of games:
 * 1: $v(S)=0$ if $|S|=1$; $v(S)$ is a random integer between $1$ and $100|S|$ if $2 \leq |S| < n$; $v(N)$ is a random integer between $100(n-2)$ and $100n$
 * 2: $v(S)=0$ if $|S|=1$; $v(S)$ is a random integer between $1$ and $50n$ otherwise
 * 3: $v(S)=0$ if $|S| < n-2$; $v(S)=1$ with probability $0.9$ if $n-2 \leq |S| \leq n-1$; $v(N)=1$
 * 4: $v(S)=0$ if $|S|=1$; $v(S)$ is a random integer between $1$ and $n$ otherwise
 * 5: weighted voting game for $n \geq 7$ with weights $w_j=\lfloor \frac{n-3}{2} \rfloor$ for $j=1,\dots,5$ and $w_j=1$ for $j=6, \dots, n$, and with quota $q=4w_1+n-4$. Then $v(S)=1$ if $\sum_{j \in S} w_j \geq q$ and $v(S)=0$ otherwise.

### Outputs

results.txt: containing variables
```
seed, time, iter, piv, sr (if applicable), x
```
separated by linebreaks
 * seed: seed provided in "input.txt"
 * time: computation time needed in seconds
 * iter: number of iterations needed
 * piv: number of pivots needed
 * sr: number of LPs solved in the subroutines (only for BNF, SP, PD and DP)
 * x: the nucleolus of the game (in n lines)

## Verifying the nucleolus

### Algorithms

For each algorithm covered there is an executable and the source as well.

Verifying folder:
 * Kohlberg: based on the original Kohlberg criterion based on [Kohlberg (1971) - On the nucleolus of a characteristic function game](https://epubs.siam.org/doi/pdf/10.1137/0120009)
 * SKA: simplified Kohlberg algorithm from Benedek et al. (2018) - Finding and verifying the nucleolus of cooperative games (forthcoming)
 * SKAcr: simplified Kohlberg algorithm with compact representation from Benedek et al. (2018) - Finding and verifying the nucleolus of cooperative games (forthcoming)

### Executables

All executables are obtained by compiling the corresponding source (see below) using Microsoft Visual Studio Enterprise 2015 (version 14.0.25431.01 Update 3) on OS Microsoft Windows 10 Enterprise (version 10.0.17134 Build 17134). NOTE THAT ALL EXECUTABLES ARE USING CPLEX (version 12.7.0)

Verifying:
 * Kohlberg.exe
 * SKA.exe
 * SKAcr.exe

### Sources

Verifying: each of the codes below requires common.cpp, common.h, gen_game.cpp and gen_game.h
 * Kohlberg.cpp, Kohlberg.h
 * SKA.cpp, SKA.h
 * SKAcr.cpp, SKAcr.cpp

### Inputs

input.txt: containing variables
```
n, type, seed, disp, memo
```
separated by linebreaks
 * n: integer, the number of players
 * type: integer between 0 and 7; type=0 reads the game from "v.txt" (see below), type>0 generates certain types of games (see below)
 * seed: integer, seed for randomization when generating certain types of games
 * disp: boolean, 1 to display information while running (default 0)
 * memo: boolean, 1 to switch to memory-saving implementation (default 0); as a rough guide, using 16 Gb memory only memo=1 works with ~n>28 for most games

v.txt: being read only if type=0
 * containing coalitional values of a game separated by linebreaks
 * for an n-player game it should contain 2^n-1 lines
 * the coalitional value in the i-th row (for i between 1 and 2^n-1) corresponds to the coalition with n-dim. characteristic vector as the binary number converted from the decimal number i
 * for example if n=5, then the file should contain 31 lines
 * the first line should contain the value of player 1
 * the 7-th line should contain the value of coalition containing players 1, 2 and 3
 * the 19-th line should contain the value of coalition containing players 1, 2 and 5

types of games:
 * 1: $v(S)=0$ if $|S|=1$; $v(S)$ is a random integer between $1$ and $100|S|$ if $2 \leq |S| < n$; $v(N)$ is a random integer between $100(n-2)$ and $100n$
 * 2: $v(S)=0$ if $|S|=1$; $v(S)$ is a random integer between $1$ and $50n$ otherwise
 * 3: $v(S)=0$ if $|S| < n-2$; $v(S)=1$ with probability $0.9$ if $n-2 \leq |S| \leq n-1$; $v(N)=1$
 * 4: $v(S)=0$ if $|S|=1$; $v(S)$ is a random integer between $1$ and $n$ otherwise
 * 5: weighted voting game for $n \geq 7$ with weights $w_j=\lfloor \frac{n-3}{2} \rfloor$ for $j=1,\dots,5$ and $w_j=1$ for $j=6, \dots, n$, and with quota $q=4w_1+n-4$. Then $v(S)=1$ if $\sum_{j \in S} w_j \geq q$ and $v(S)=0$ otherwise.

sol.txt:
 * the solution to test (in n lines), coorinates separated by linebreaks

### Outputs

results.txt: containing variables
```
seed, dec, time, iter, sr, saved (if applicable), x
```
separated by linebreaks
 * seed: seed provided in "input.txt"
 * dec: boolean, 1 if the tested solution is the nucleolus
 * time: computation time needed in seconds
 * iter: number of iterations needed
 * sr: number of LPs solved in the subroutines
 * saved: number of settled coalitions not needed to store (only for SKAcr)
 * x: tested solution provided in "sol.txt" (in n lines)

## Testgames

An extensive amount of test instances are provided as well.

v$i$n$k$j$l$.txt:
 * corresponds to the $l$-th $i$-player type $k$ game
 * for $1 \leq i \leq 4$ and $k \in \{5,10,15,20,25,26,27,28\}$, $l$ goes from $1$ to $50$
 * for $1 \leq i \leq 4$ and $29 \leq k \leq 30$, $l$ goes from $1$ to $10$
 * for $i=5$ and $k \in \{7,10,15,20,25,26,27,28\}$, $l=1$ (the game is deterministic)
 * each file consists of $4n+5$ lines
 * the first line contains the seed, in the next $4(n+1)$ four solutions are given
 * the four solutions are a random imputation, an element of the least core, an element of the $leastest core$ (solving the first two LPs in the LP sequence) and the nucleolus
 * each solution (in n lines) is followed by a boolean, taking the value 1 if and only if the corresponding solution is the nucleolus

## Examples

To find the nucleolus of the following 3-player game
```
$v(\{1\})=1, v(\{2\})=3, v(\{3\})=5, v(\{1,2\})=6, v(\{1,3\})=7, v(\{2,3\})=8, v(N)=12$
```
the file "v.txt" should contain
```
1
2
6
5
7
8
12
```
While the "input.txt" should consist of
```
3
0
0
1
0
0
```
if we would like to get information about the progress of the chosen algorithm (the seed can be arbitrary if type=0). Then we can run "BNF.exe" and obtain the solution
```
2.75
3.75
5.5
```
If we save this in "sol.txt" we can verify the solution using one of the verification algorithms.

However if we were to find the nucleolus of a 15-player type 2 game using seed 480453, without displaying information but with limiting memory usage and without linear speed-up we should change input.txt to
```
15
2
480453
0
1
1
```

## Authors

* **Marton Benedek** - *Initial work* - [Nucleolus](https://github.com/blrzsvrzs/nucleolus)

## License

This project is licensed under the GNU General Public License - see the [LICENSE](LICENSE) file for details
