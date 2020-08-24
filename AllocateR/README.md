This folder contains Allocation Files that take advantage of the PrjOptiAlloc Project.

1. First, open up three files in alloc_discrete_fun, feasible, feasible_threshold as well as optimal, all three files
2. Go to L63 for each of the file, update CSV name.
3. Go to title, change title
4. Update bin width or other parameters, and also note in title
5. go back to: \PrjOptiSNW\AllocateR\fs_vig_render.Rmd

## The state-space of the Underlying Model

We solve for the value of 244 checks, each at $100, at these state-space points. We computed the marginal effects on consumption and value of 1st to 244th check for each element of the state-space:

- 2 marital groups: Marital Status 0 or 1
- 5 kids groups: Children 0 to 4
- 2 education groups: college or not
- 83 age groups: Age 18 to 100
- 65 savings groups
- 266 household head productivity shock.
- 5 spousal productivity shocks if married

Overall state-space: 86,104,200
  + married state-space: 5*2*83*65*266*5 = 71,753,500
  + single state-space (no spousal shock): 5*2*83*65*266*1 = 14,350,700

## 2020 Unemployment Shock COVID

There is a covid unemployment shock in 2020, that lowers (proportional) income realization.

So we solve the 86,104,200 twice once with the unemployment shock, once without

- State-space points total with COVID shock: 172,208,400

Since we have 244 checks, we solve for the consumption increase and value increase from $100 increment of the check 244 times, + 1 time for without check.

So overall, we have 42.2 billion MPCs: (244+1)*172208400=42,191,058,000

## How Many Cells We Solved for from 2019 Perspective for the COVID planner

The planner does not observe the full state-space so we intergrate from 2019 perspective given 2019 to 2020 COVID transition probability (conditional on state-space) and kids and shock transitions planner value.

- 5 kids groups: Children 0 to 4
- 2 marital groups: Marital Status 0 or 1
- 83 age groups: Age 18 to 100
- 509 Income Groups:
    + 476 bins below max actual phaseout: solved at $500 intervals between $0 and $238,000
    + 33 bins after max actual phaseout: solved at $5000 interval after $238,000, where the 33rd final bin is between $401,130 and Maximum
    + The income group is composed of individuals of with different 2019 productivity shocks and savings levels.

Together, these are:

- 5*2*83*509 = 422,470 groups/bins

We have (244+1) marginal average consumption gain, and value gain for each of the groups, so from 2019 planner perspective, we have 103.5 million MPCs

- 103,505,150 = 422470*(244+1)

## How Many Cells do we consider for the Allocation Problem

For the optimal allocation problem, we have several versions, how many individual types are there in each version?

- for FEASIBLE allocation, there are 970=5*2*97 types/cells of households:
    + 5 children groups
    + 2 spousal groups
    + 97 income bins: the allocation planner sees approximately $2500 income bins between $0 and $238,000, and 1 bin after $238,000. There are 97 bins
- for OPTIMAL G4 (4 age groups 18 to 64) allocation, there are 3880=5*2*97*4 types/cells of households:
    + 5 children groups
    + 2 spousal groups
    + 4 age groups
    + 97 income bins
- for OPTIMAL G47 (47 age groups) allocation, there are 45590=5*2*97*47 types/cells of households:
    + 5 children groups
    + 2 spousal groups
    + 47 age groups
    + 97 income bins
