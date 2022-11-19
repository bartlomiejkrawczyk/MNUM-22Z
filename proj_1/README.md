# Metody Numeryczne - Projekt I
```s
student: Bartłomiej Krawczyk
indeks: 310774
```

# Dane

### A):
$$
\begin{equation}
  a_{ij} = \begin{cases}
    -10, & \text{dla $j = i$} \\
    3, & \text{dla $j = i - 1$ lub $j = i + 1$} \\
    0, & \text{dla pozostałych}
  \end{cases}
\end{equation}
$$

$$
\begin{equation}
  b_{ij} = 2.5 - 0.5i
\end{equation}
$$

**Funkcja generująca macierz A oraz wektor b:**

```matlab
function [A,b] = prepareParametersA(n)
    A = zeros(n, n);
    b = zeros(n, 1);

    for i = 1 : n
        A(i, i) = -10;
        b(i, 1) = 2.5 - 0.5 * i;
    end
    
    for i = 2 : n
        A(i, i - 1) = 3;
        A(i - 1, i) = 3;
    end
end
```

Przykładowa macierz $A$ dla $n = 5$:
```s
   -10     3     0     0     0
     3   -10     3     0     0
     0     3   -10     3     0
     0     0     3   -10     3
     0     0     0     3   -10
```

Przykładowy wektor $b$ dla $n = 5$:
```s
    2.0000
    1.5000
    1.0000
    0.5000
         0
```

### B):
$$
\begin{equation}
  a_{ij} = \begin{cases}
    4n^2 + (2i + 3) n, & \text{dla $j = i$} \\
    2 (i + j) + 1, & \text{dla $j \neq i $} \\
  \end{cases}
\end{equation}
$$

$$
\begin{equation}
  b_{ij} = 2.5 + 0.6i
\end{equation}
$$

**Funkcja generująca macierz A oraz wektor b:**

```matlab
function [A,b] = prepareParametersB(n)
    A = zeros(n, n);
    b = zeros(n, 1);

    for i = 1 : n
        for j = 1 : n
            A(i, j) = 2 * (i + j) + 1;
        end
    end
    
    for i = 1 : n
        A(i, i) = 4 * n ^ 2 + (2 * i + 3) * n;
        b(i, 1) = 2.5 + 0.6 * i;
    end
end
```

Przykładowa macierz $A$ dla $n = 5$:
```s
   125     7     9    11    13
     7   135    11    13    15
     9    11   145    15    17
    11    13    15   155    19
    13    15    17    19   165
```

Przykładowy wektor $b$ dla $n = 5$:
```s
    3.1000
    3.7000
    4.3000
    4.9000
    5.5000
```

# Zadanie 1

## Treść

Napisać uniwersalną procedurę w Matlabie o odpowiednich parametrach wejścia i wyjścia (solwer), rozwiązującą układ $n$ równań liniowych $Ax = b$, wykorzystując podaną metodę.

> Nie sprawdzać w procedurze, czy dana macierz $A$ spełnia wymagania stosowalności metody.

Obliczyć błąd rozwiązania $\varepsilon = ∥A\tilde{x} − b∥_2$ (skorzystać z funkcji `norm` Matlaba).

Proszę zastosować następnie swoją procedurę w programie do rozwiązania obydwu (jeśli można) lub jednego z układów równań dla podanych niżej macierzy $A$ i wektorów $b$, przyjmując: $n = 5, 10, 25, 50, 100, 200$.

Metoda: faktoryzacji $LDL^T$

Proszę wykonać wykres (wykresy) zależności błędu $\varepsilon$ od liczby równań $n$.

## Rozwiązanie



$$
\begin{align}
y = a
\end{align}
$$

```matlab
function x = solve_linar_equation_with_triangular_matrix_upper(A, b)

    [n, ~] = size(A);
    x = zeros(n, 1);
    
    for k = n : -1 : 1
        x(k, 1) = b(k, 1);
        
        for j = k + 1 : n
            x(k, 1) = x(k, 1) - A(k, j) * x(j, 1);
        end

        x(k, 1) = x(k, 1) / A(k, k);
    end
end
```


# Zadanie 2

## Treść

Napisać uniwersalną procedurę w Matlabie o odpowiednich parametrach wejścia i wyjścia, rozwiązującą układ $n$ równań liniowych $Ax = b$, wykorzystując metodę **iteracyjną Jacobiego**.

> Nie sprawdzać w procedurze, czy dana macierz $A$ spełnia wymagania stosowalności metody.

Jej parametry wejściowe powinny zawierać m.in. wartość graniczną $\delta$ błędu między kolejnymi przybliżeniami rozwiązania, liczonego jako **norma euklidesowa z ich różnicy** (skorzystać z funkcji `norm` Matlaba). Przyjąć jako kryterium stopu warunek: δ = 10−8 ≜ 1e − 8.

$$
\delta = 10^{-8} \triangleq 1e-8
$$

Proszę zastosować tę procedurę do rozwiązania właściwego układu równań spośród przedstawionych poniżej dla $n = 5, 10, 25, 50, 100, 200$.

Proszę sprawdzić dokładność rozwiązania licząc także błąd $\varepsilon$ i dla każdego układu równań wykonać rysunek zależności tego błędu od liczby równań $n$. Jeśli był rozwiązywany ten sam układ równań, co w p. 1, proszę porównać czasy obliczeń dla różnych algorytmów i wymiarów zadań.

## Rozwiązanie

# Zadanie 3

## Treść

Dla podanych w tabeli danych pomiarowych (próbek) metodą najmniejszych kwadratów należy wyznaczyć funkcję wielomianową y = f(x) (tzn. wektor współczynników) najlepiej aproksymującą te dane.

| $x_i$ | $y_i$   |
|-------|---------|
| -10   | -42.417 |
| -8    | -23.440 |
| -6    | -11.160 |
| -4    | -4.128  |
| -2    | -0.725  |
| 0     | 0.942   |
| 2     | -2.069  |
| 4     | -3.908  |
| 6     | -4.705  |
| 8     | -5.438  |
| 10    | -3.578  |

Proszę przetestować wielomiany stopni: $3, 5, 7, 9, 10$. Kod aproksymujący powinien być uniwersalną procedurą w Matlabie o odpowiednich parametrach wejścia i wyjścia.

W sprawozdaniu proszę przedstawić na rysunku otrzymaną funkcję na tle danych (funkcję aproksymującą proszę próbkować przynajmniej 10 razy częściej niż dane).

Do rozwiązania zadania najmniejszych kwadratów proszę wykorzystać najpierw **układ równań normalnych**, a potem **rozkład SVD**.

Do rozwiązywania układu równań i dekompozycji użyć solwerów Matlaba.. Porównać efektywność obydwu podejść.

Do liczenia wartości wielomianu użyć funkcji `polyval`.

Proszę obliczyć błąd aproksymacji w dwóch normach: euklidesowej oraz maksimum (nieskończoność). W obydwu przypadkach skorzystać z funkcji `norm` Matlaba.

## Rozwiązanie

