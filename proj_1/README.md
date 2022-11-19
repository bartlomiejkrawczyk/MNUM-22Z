# Metody Numeryczne - Projekt I
```s
student: Bartłomiej Krawczyk
indeks: 310774
```


# Zadanie 1

Napisać uniwersalną procedurę w Matlabie o odpowiednich parametrach wejścia i wyjścia (solwer), rozwiązującą układ $n$ równań liniowych $Ax = b$, wykorzystując podaną metodę. Nie sprawdzać w procedurze, czy dana macierz $A$ spełnia wymagania stosowalności metody. Obliczyć błąd rozwiązania $\varepsilon = ∥A\tilde{x} − b∥_2$ (skorzystać z funkcji `norm` Matlaba).

Proszę zastosować następnie swoją procedurę w programie do rozwiązania obydwu (jeśli można) lub jednego z układów równań dla podanych niżej macierzy $A$ i wektorów $b$, przyjmując: $$n = 5, 10, 25, 50, 100, 200$$

Metoda: faktoryzacji $LDL^T$

Proszę wykonać wykres (wykresy) zależności błędu $\varepsilon$ od liczby równań $n$.




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


# Zadanie 3



