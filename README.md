## FOREST_FIRE_SIMULATION

Approximated rate of expansion and pattern of a pine tree forest fire at any coordinates.

It is based on the cellular automata type code with a few changes in process.

Takes into consideration: weather (humidity level, wind speed and angle, temperature), elevation difference between cells and average humidity level of a pine tree. Each based on the inputted coordinates.

*Inputs:* Latitude coordinates, Longitude coordinates, Grid size, Time that passes between each grid display.

*Outputs:* Introduction, (repeated:) Time (seconds, minutes, hours, days) passed, Grid; a visualization of the fire spread as a 2D list of lists.

[Formulas source](https://www.mdpi.com/1999-4907/14/7/1371#B62-forests-14-01371)
