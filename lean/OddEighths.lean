import Mathlib.Tactic.Ring
import Mathlib.Algebra.Ring.Basic
import Mathlib.Data.Int.GCD
import Mathlib.Data.List.Basic

/-!
# La espina aritmetica de la Parte VII, y solo la espina

Copyright (c) 2026 Carles Marin. Author: Carles Marin (with Claude, Anthropic, as assistant).

## Que se formaliza aqui, y por que solo esto

La cadena de la Parte VII se come un coeficiente que el articulo NO calcula: el aporte del sector
gauge, leido de la ec. (68) de Komori-Maru. Ese numero esta en cuestion. Lean **no puede
adjudicarlo** -- es un input fisico, no un paso logico -- asi que formalizar el reparto
`(3/2, 1/2)` seria justamente adjudicar la ecuacion de otro, y eso no se hace aqui.

Lo que si se puede blindar, y es lo que hay debajo, es la parte que **no depende de la semilla**:

* `mod6_law`  : `8D = 2*A4 + 3 (mod 6)` para todo contenido -- el Teorema 2 del articulo.
* `twoW_odd`  : `2W` es impar para todo contenido.
* `oddEighths`: el Teorema 1, en su forma **condicional** -- si el aporte gauge es impar,
                entonces `8D` es impar.

Las tres valen con las DOS semillas, y eso es el hallazgo que las hace formalizables ahora: el
punto base gauge da el mismo resto `9 = 3 (mod 6)` en los dos casos, de modo que la ley mod 6 no
distingue las ramas. El Teorema 1 si las distingue, y por eso se enuncia con su hipotesis fuera.

## La eleccion de coordenadas, que no es cosmetica

Se trabaja en `(2*A4, 8D, 2W)` y no en `(A4, 8D, 2W)`. Con la semilla candidata `A4` es
SEMI-ENTERO, asi que `A4` no es una coordenada entera y `2*A4` si. Y la combinacion en la que la
prueba sale sin mencionar ninguna paridad, `8D - 2*A4`, es exactamente la que Lean quiere. Que la
prueba correcta y la coordenada comoda coincidan suele ser senal de que la estructura es la buena.
-/

namespace PartVII

/-- Los ocho multipletes del bulk, en coordenadas `(2*A4, 8D, 2W)`.
    Leidos de `outputs/gauge_ghost_seed.txt`, control C5. `A4` es entero para la materia, de modo
    que la primera coordenada es su doble. -/
def matter : List (ℤ × ℤ × ℤ) :=
  [ (0,   -6,   2),      -- 7(+,+)
    (2,    8,  -2),      -- 7(+,-)
    (34,  16,   6),      -- 28(+,+)
    (8,    2,  -6),      -- 28(+,-)
    (36,   0,  12),      -- 48(+,+)
    (16,  28, -12),      -- 48(+,-)
    (136, 10,  18),      -- 84(+,+)
    (218, 80, -18) ]     -- 84(+,-)

/-- El sector gauge, en las mismas coordenadas. Las dos semillas de la bifurcacion. -/
def gaugePublished : ℤ × ℤ × ℤ := (-36, -27, -3)   -- la ec. (68) tal y como se imprime
def gaugeCandidate : ℤ × ℤ × ℤ := (-27, -18, -3)   -- con la resta de fantasmas

/-- La combinacion que lo decide todo: `8D - 2*A4`. -/
def key (v : ℤ × ℤ × ℤ) : ℤ := v.2.1 - v.1

/-- La tercera coordenada, `2W`. -/
def w2 (v : ℤ × ℤ × ℤ) : ℤ := v.2.2

-- ## Ladrillo 1 -- la materia, comprobada uno a uno

/-- Todo multiplete de materia cumple `8D - 2*A4 = 0 (mod 6)`. Comprobacion finita. -/
theorem matter_key_mod6 : ∀ v ∈ matter, (6 : ℤ) ∣ key v := by decide

/-- Todo multiplete de materia aporta `2W` PAR. -/
theorem matter_w2_even : ∀ v ∈ matter, (2 : ℤ) ∣ w2 v := by decide

-- ## Ladrillo 3 -- el punto base, y EL MISMO 9 en las dos ramas

/-- Con la semilla publicada el punto base da `9`. -/
theorem key_published : key gaugePublished = 9 := by decide

/-- Con la semilla candidata da **el mismo** `9`. Este es el hecho que hace que el Teorema 2
    sobreviva a la bifurcacion mientras el Teorema 1 no. -/
theorem key_candidate : key gaugeCandidate = 9 := by decide

/-- Dicho de una vez: las dos semillas son indistinguibles para `8D - 2*A4`. -/
theorem key_seed_independent : key gaugePublished = key gaugeCandidate := by decide

/-- Y las dos aportan `2W` impar. -/
theorem gauge_w2_odd : ¬ (2 : ℤ) ∣ w2 gaugePublished ∧ ¬ (2 : ℤ) ∣ w2 gaugeCandidate := by decide

-- ## Ladrillo 2 -- la linealidad, que es donde entra el contenido arbitrario

/-- Un contenido es una lista de multiplicidades no negativas contra los generadores.
    `combo ns vs` es la suma ponderada de una coordenada. -/
def combo (ns : List ℕ) (xs : List ℤ) : ℤ :=
  (List.zipWith (fun (k : ℕ) (x : ℤ) => (k : ℤ) * x) ns xs).sum

/-- Si `d` divide a cada generador, divide a cualquier combinacion con multiplicidades naturales.
    Es el paso que lleva de los ocho casos comprobados a TODO contenido. -/
theorem dvd_combo {d : ℤ} (ns : List ℕ) (xs : List ℤ) (h : ∀ x ∈ xs, d ∣ x) :
    d ∣ combo ns xs := by
  unfold combo
  induction ns generalizing xs with
  | nil => simp
  | cons k ks ih =>
    cases xs with
    | nil => simp
    | cons x xs =>
      -- se construye el testigo a mano en vez de invocar `dvd_add`: en este Mathlib ese nombre
      -- no resuelve con imports especificos, y el testigo explicito no necesita ningun lema.
      simp only [List.zipWith_cons_cons, List.sum_cons]
      obtain ⟨c, hc⟩ := h x (by simp)
      obtain ⟨e, he⟩ := ih xs (fun y hy => h y (by simp [hy]))
      exact ⟨(k : ℤ) * c + e, by rw [hc, he]; ring⟩

-- ## Los tres enunciados

/-- **Teorema 2, y es independiente de la semilla.** Para cualquier punto base gauge cuyo
    `8D - 2*A4` valga `9` -- lo cual cubre las DOS semillas por `key_published` y
    `key_candidate` -- y cualquier contenido con multiplicidades naturales,
    `8D - 2*A4 = 3 (mod 6)`, es decir `8D = 2*A4 + 3 (mod 6)`. -/
theorem mod6_law (g : ℤ × ℤ × ℤ) (hg : key g = 9) (ns : List ℕ) :
    (key g + combo ns (matter.map key)) % 6 = 3 := by
  have hc : (6 : ℤ) ∣ combo ns (matter.map key) := by
    refine dvd_combo ns _ ?_
    intro x hx
    obtain ⟨v, hv, rfl⟩ := List.mem_map.mp hx
    exact matter_key_mod6 v hv
  obtain ⟨t, ht⟩ := hc
  rw [hg, ht]
  omega

/-- **El teorema de `2W`, tambien independiente de la semilla.** El sector gauge aporta `-3`,
    impar, con las dos; la materia aporta par; luego `2W` es impar para todo contenido. -/
theorem twoW_odd (g : ℤ × ℤ × ℤ) (hg : w2 g = -3) (ns : List ℕ) :
    ¬ (2 : ℤ) ∣ (w2 g + combo ns (matter.map w2)) := by
  have hc : (2 : ℤ) ∣ combo ns (matter.map w2) := by
    refine dvd_combo ns _ ?_
    intro x hx
    obtain ⟨v, hv, rfl⟩ := List.mem_map.mp hx
    exact matter_w2_even v hv
  obtain ⟨t, ht⟩ := hc
  rw [hg, ht]
  omega

/-- **Teorema 1, en su forma condicional.** La materia aporta par a `8D`; asi que SI el aporte
    gauge es impar, `8D` es impar para todo contenido -- y entonces `8D` no puede ser cero.

    La hipotesis va fuera a proposito. Es el unico eslabon que este articulo no calcula: se lee
    de la ec. (68) de Komori-Maru. Si la semilla se mueve, el teorema SIGUE SIENDO CIERTO y su
    hipotesis deja de cumplirse, que es una situacion distinta de un teorema que falla. -/
theorem oddEighths (dGauge : ℤ) (hodd : ¬ (2 : ℤ) ∣ dGauge) (ns : List ℕ)
    (hmatter : (2 : ℤ) ∣ combo ns (matter.map (fun v => v.2.1))) :
    ¬ (2 : ℤ) ∣ (dGauge + combo ns (matter.map (fun v => v.2.1))) := by
  obtain ⟨t, ht⟩ := hmatter
  rw [ht]
  omega

/-- La materia aporta par a `8D`, que es la hipotesis que `oddEighths` necesita y que aqui se
    comprueba en vez de suponerse. -/
theorem matter_d_even : ∀ v ∈ matter, (2 : ℤ) ∣ v.2.1 := by decide

/-- **Corolario para la semilla publicada.** Su `8D` gauge es `-27`, impar, luego la hipotesis
    del Teorema 1 se cumple y `8D` es impar para todo contenido. -/
theorem oddEighths_published (ns : List ℕ) :
    ¬ (2 : ℤ) ∣ ((-27 : ℤ) + combo ns (matter.map (fun v => v.2.1))) := by
  refine oddEighths (-27) (by decide) ns ?_
  refine dvd_combo ns _ ?_
  intro x hx
  obtain ⟨v, hv, rfl⟩ := List.mem_map.mp hx
  exact matter_d_even v hv

/-- **Y para la candidata, la hipotesis NO se cumple**: su `8D` gauge es `-18`, par. El teorema
    no falla; su premisa no se satisface. -/
theorem candidate_hypothesis_fails : (2 : ℤ) ∣ (gaugeCandidate.2.1) := by decide

-- El certificado, y no es decorativo: `compila` no es lo mismo que `sin sorry y sin axiomas
-- extra`.  Lo esperado son a lo sumo los tres estandar de Lean -- propext, Classical.choice,
-- Quot.sound -- y ningun `sorryAx`.  Si aparece `sorryAx`, el ladrillo no vale.
#print axioms mod6_law
#print axioms twoW_odd
#print axioms oddEighths
#print axioms oddEighths_published
#print axioms key_seed_independent

end PartVII
