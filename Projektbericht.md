# Projektbericht
Multiplikation dünnbesetzter Matrizen in CSR Format

## Aufgabenstellung und Entwicklung

Die Aufgabe dieses Projekts ist es, eine effiziente Multiplikation dünnbesetzter Matrizen umzusetzen.
Hierfür muss CSR verwendet, der Multiplikationsalgorithmus entsprechend angepasst und in ein Rahmenprogramm eingefügt werden.
V1 wurde als Testimplementierung mit einer Konvertierung von und zu CSR und der sequentiellen Multiplikation von vollbesetzten Matrizen entwickelt.
Der Gustavson-Algorithmus ist die Basis für die Multiplikation ((Vergleichs-) Implementierungen ab V2) und trägt maßgeblich zur Effizienz der Multiplikation durch Ignorierung der Nullwerte bei.
Speichereffizienz wird durch die Abschätzung der Werte der Ergebnismatrix erreicht.
In V3 (SSE) und V4 (AVX) wird SIMD umgesetzt. In V0 wird zur Laufzeit entschieden, ob Multithreading verwendet wird und wie viele Threads dem Prozess alloziert werden sollen. Alternativ wird V6 aufgerufen. 

## Implementierungen
V0: Gustavson-Algorithmus, Multithreading, Strategy Pattern
→ Hauptimplementierung
V1: Konvertierung (CSR ↔ dichtes Format), sequentielle (vollbesetzte) Multiplikation
V2: Gustavson-Algorithmus →  Ausnutzung des CSR-Formats
V3: Gustavson-Algorithmus, SIMD (SSE)
V4: Gustavson-Algorithmus, SIMD (AVX)
V5: Gustavson-Algorithmus, Größenabschätzung

## Benchmarking
Getestet wurde auf einer geeigneten Linux-Maschine geringer Auslastung durch andere Prozesse, mit zufällig durch Seed generierten Matrizen unterschiedlicher Größe und Dichte, wobei kleinere Berechnungen mehrmals ausgeführt wurden, um aussagekräftige Benchmarks zu erhalten.
- Intel(R) Core(TM) i5-4460  CPU @ 3.20GHz, 4-Kerne, 64-bit, SSE 4.2, AVX2
- 32 GiB RAM
- Ubuntu 23.10, Linux 6.5.0-42-generic

## Einzelleistungen

Berk Erdemoglu verknüpfte das Programm zu einem Ganzen, indem er das Rahmenprogramm mit E/A und die Main-Methode mit Argument Parsing erstellte.
Darüber hinaus entwickelte er das Multithreading aus der Hauptimplementierung und die SSE-Implementierung (V3). Schließlich dokumentierte Berk das Projekt mit docstrings.
Mit seinen fortgeschrittenen C-Kenntnissen half er dem Team in unterschiedlichen Problemstellungen aus und vereinigte das Programm.

Erik Eremenko recherchierte gemeinsam mit Alessandro mögliche Optimierungsstrategien. Er implementierte die erste Matrixmultiplikation (V1) als Vergleichsimplementierung zum Testen und entwickelte ein Generatorprogramm zur Erstellung von Testmatrizen, was für realistische Tests genutzt wurde. Des Weiteren erstellte er diverse Hilfsmethoden und führte gemeinsam mit Berk das Benchmarking durch. Er entwickelte in V6 die Bestimmung der Anzahl der zu verwendenden Threads und erstellte die Präsentation.

Alessandro Kronast-Reichert trug der Recherche maßgeblich mit dem Vorschlag zur Nutzung des Gustavson-Algorithmus bei und implementierte ihn für die Multiplikationen (V0/V2/V3). Er stellte die Wichtigkeit der Abschätzung der Größe der Ergebnismatrix fest und implementierte sie, was die Performanz erheblich verbesserte. Außerdem trug er dem Testen der Implementierungen mit einer Methode bei, welche die Gleichheit von CSR-Matrizen mit unterschiedlichen Column-Value-Permutationen prüfen kann.

---

Anmerkung: Dieser Projektbericht wurde für die Veröffentlichung bearbeitet und Details der Aufgabenstellung zur Modulprüfung entfernt.
