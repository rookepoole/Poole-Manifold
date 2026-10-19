(* ::Package:: *)
(* Poole Manifold | PM0 base cellular automaton | Wolfram Community research companion *)
(* Rooke Alan Poole: originating research. ChatGPT-assisted translation and test harness. *)
(* 18 October 2026 | Wolfram Language | Exact finite integer dynamics *)
(* All calculations use finite periodic cubic arrays; state x(t) is Markovian. *)
(* This code DOES NOT implement prime-sharpening, quantum dynamics, or physical models. *)

ClearAll[
  pmRule, pmAblatedRule, pmOffsets3, pmOffsets2, pmCount3,
  pmNativeStep, pmReferenceStep, pmRectMask, pmSheet, pmRectState,
  pm2DCounts, pmPlanarPredict, pmRectangleTrial, pmRandomMaskChecks,
  pmRandomFull3DChecks, pmRunChecks, pmFigure3D, pmFigure2D,
  pmExample, pmNext, pmBorn, pmZ, pmSweep, pmEvolution
];

(* Section 1 | Native Wolfram rule and independent cyclic reference *)
pmRule = <|
  "Dimension" -> 3,
  "Neighborhood" -> "Moore",
  "GrowthSurvivalCases" -> {{5, 6, 7}, {5, 6, 7, 8, 9}}
|>;

pmAblatedRule = <|
  "Dimension" -> 3,
  "Neighborhood" -> "Moore",
  "GrowthSurvivalCases" -> {{5, 6, 7}, {5, 6, 7}}
|>;

pmOffsets3 = DeleteCases[Tuples[{-1, 0, 1}, 3], {0, 0, 0}];
pmOffsets2 = DeleteCases[Tuples[{-1, 0, 1}, 2], {0, 0}];

(* RotateLeft with a 3-component shift wraps periodically along all axes. *)
pmCount3[x_] := Total[RotateLeft[x, #] & /@ pmOffsets3];
pmNativeStep[x_] := Last[CellularAutomaton[pmRule, x, 1]];
pmReferenceStep[x_] := Module[{n = pmCount3[x]},
  UnitStep[n - 5] UnitStep[7 + 2 x - n]
];

(* Section 2 | Sheet, defect geometry, exact q/r predictor *)
(* Inputs d are square 0/1 masks; 1 marks a missing site in the sheet. *)
pmSheet[d_] := Module[{l = Length[d], z, x},
  z = Quotient[l + 1, 2];
  x = ConstantArray[0, {l, l, l}];
  x[[z]] = 1 - d;
  x
];

(* Place a rectangle inside a lattice with a two-voxel nonwrapping margin. *)
pmRectMask[l_Integer, a_Integer, b_Integer] := Module[{d},
  If[!(l >= 7 && 2 <= a <= l - 7 && 2 <= b <= l - 8),
    Return[$Failed]
  ];
  d = ConstantArray[0, {l, l}];
  d[[5 ;; 4 + a, 6 ;; 5 + b]] = 1;
  d
];

pmRectState[a_Integer, b_Integer, l_Integer : 25] :=
  pmSheet[pmRectMask[l, a, b]];

pm2DCounts[d_] := Module[{q},
  q = Total[RotateLeft[d, #] & /@ pmOffsets2];
  {q, q + d}
];

(* Exact one-step prediction, expressed purely in the planar defect counts. *)
pmPlanarPredict[d_] := Module[{l = Length[d], z, q, r, core, halo, y},
  z = Quotient[l + 1, 2];
  {q, r} = pm2DCounts[d];
  core = (1 - d) UnitStep[3 - q] + d UnitStep[q - 1] UnitStep[3 - q];
  halo = UnitStep[r - 2] UnitStep[4 - r];
  y = ConstantArray[0, {l, l, l}];
  y[[z]] = core;
  y[[z - 1]] = halo;
  y[[z + 1]] = halo;
  y
];

(* Section 3 | First-step rectangular-hole counts and typed birth channels *)
pmRectangleTrial[a_Integer, b_Integer, l_Integer : 25] :=
  Module[{d, x, y, n, born, neighborsAtBirth, actualChannels},
    d = pmRectMask[l, a, b];
    If[d === $Failed, Return[$Failed]];
    x = pmSheet[d];
    y = pmNativeStep[x];
    n = pmCount3[x];
    born = (1 - x) y;
    neighborsAtBirth = Pick[Flatten[n], Flatten[born], 1];
    actualChannels = Table[Count[neighborsAtBirth, k], {k, 5, 7}];
    <|
      "a" -> a, "b" -> b,
      "MeasuredBirths" -> Length[neighborsAtBirth],
      "FormulaBirths" -> 4 a + 4 b + 12,
      "MeasuredDeaths" -> Total[Flatten[x (1 - y)]],
      "MeasuredN5N6N7" -> actualChannels,
      "FormulaN5N6N7" -> {12, 4 a + 4 b - 16, 16},
      "FullPlanarPredictionMatches" -> (y === pmPlanarPredict[d]),
      "NativeVsShiftReferenceMatches" -> (y === pmReferenceStep[x])
    |>
  ];

(* Section 4 | Deterministic checks for local methods and rule ablation *)
pmRandomMaskChecks[count_Integer : 12, l_Integer : 17, seed_Integer : 20261018] :=
  BlockRandom[
    SeedRandom[seed];
    Table[
      Module[{d = RandomInteger[1, {l, l}]},
        pmNativeStep[pmSheet[d]] === pmPlanarPredict[d]
      ],
      {count}
    ]
  ];

pmRandomFull3DChecks[count_Integer : 5, l_Integer : 9, seed_Integer : 20261018] :=
  BlockRandom[
    SeedRandom[seed];
    Table[
      Module[{x = RandomInteger[1, {l, l, l}]},
        pmNativeStep[x] === pmReferenceStep[x]
      ],
      {count}
    ]
  ];

pmRunChecks[] := Module[
  {l = 25, z = 13, plane, pore, ablated, reports, rectPass},
  plane = pmSheet[ConstantArray[0, {l, l}]];
  pore = ConstantArray[0, {l, l}]; pore[[10, 10]] = 1;
  ablated = Last[CellularAutomaton[pmAblatedRule, plane, 1]];
  reports = pmRectangleTrial @@@ {{2, 2}, {3, 4}, {10, 10}, {12, 12}};
  rectPass = And @@ (Function[rec,
       rec["MeasuredBirths"] === rec["FormulaBirths"] &&
       rec["MeasuredDeaths"] === 0 &&
       rec["MeasuredN5N6N7"] === rec["FormulaN5N6N7"] &&
       TrueQ[rec["FullPlanarPredictionMatches"]] &&
       TrueQ[rec["NativeVsShiftReferenceMatches"]]
     ] /@ reports);
  <|
    "Moore26" -> (Length[pmOffsets3] === 26),
    "StableFullPlane" -> (pmNativeStep[plane] === plane),
    "StableSinglePore" -> (pmNativeStep[pmSheet[pore]] === pmSheet[pore]),
    "AblationErasesSheet" -> (Total[Flatten[ablated]] === 0),
    "SampleRectanglesPass" -> rectPass,
    "Full3DRandomStateParity" -> And @@ pmRandomFull3DChecks[],
    "RandomPlanarMaskSelectorParity" -> And @@ pmRandomMaskChecks[],
    "RectangularReports" -> reports
  |>
];

(* Section 5 | Notebook-native figures; execute in a Wolfram front end *)
(* Cropped slab avoids plotting the full 25^3 cube. *)
pmExample = pmRectState[3, 4];
pmNext = pmNativeStep[pmExample];
pmBorn = (1 - pmExample) pmNext;
pmZ = 13;

pmFigure3D[] := GraphicsRow[
  {
    ArrayPlot3D[Take[pmExample, {pmZ - 1, pmZ + 1}, {2, 14}, {3, 16}],
      ColorRules -> {1 -> RGBColor[.35, .49, .58], 0 -> None},
      Mesh -> False, Boxed -> False, ImageSize -> 320],
    ArrayPlot3D[Take[pmBorn, {pmZ - 1, pmZ + 1}, {2, 14}, {3, 16}],
      ColorRules -> {1 -> RGBColor[.95, .61, .27], 0 -> None},
      Mesh -> False, Boxed -> False, ImageSize -> 320],
    ArrayPlot3D[Take[pmNext, {pmZ - 1, pmZ + 1}, {2, 14}, {3, 16}],
      Mesh -> False, Boxed -> False, ImageSize -> 320]
  }, Spacings -> 5, ImageSize -> 950
];

pmFigure2D[] := GraphicsRow[
  {
    ArrayPlot[pmExample[[pmZ]], Mesh -> False, ImageSize -> 300],
    ArrayPlot[pmBorn[[pmZ]], Mesh -> False, ImageSize -> 300],
    ArrayPlot[pmBorn[[pmZ + 1]], Mesh -> False, ImageSize -> 300]
  }, Spacings -> 5, ImageSize -> 950
];

(* Section 6 | Optional extended sweeps, evolution; run explicitly in notebook *)
pmSweep[low_Integer : 2, high_Integer : 9, l_Integer : 25] :=
  pmRectangleTrial @@@ Tuples[{Range[low, high], Range[low, high]}];

pmEvolution[a_Integer : 3, b_Integer : 4, t_Integer : 4] :=
  NestList[pmNativeStep, pmRectState[a, b], t];

(* No costly sweep, plot, or export happens automatically at Get time. *)
