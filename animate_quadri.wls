#!/usr/bin/env wolframscript

imageFolder = "output";
finalSize = 720;
imageSize = finalSize/3 * 6/5 * {4, 3};
exportToPov = False;
exportToPov = True;
frameCount = 180;

Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

epsilon = 1.0*^-6;
intersect[p1_, p2_, p3_, p4_] := l1 p1 + (1 - l1) p2 /. Solve[ l1 p1 + (1 - l1) p2 == l2 p3 + (1 - l2) p4, {l1, l2}][[1]];
padZero[p2d_] := Join[p2d, {0}];

(* 
the following function checks if the quadrilateral is convex and the vertices are labeled in the clockwise or counter clockwise manner. It will return false if the convex quadrilateral is like:

1 3
2 4

The probability of this event is 1/3 of the probability that the quadrilatel is convex.
1/3 * 25/36 = 25/108 = 0.231481
*)
convexThreshold = 0.001;
isConvex[quad_] := Length[NSolve[ {
        l1 quad[[1]] + (1 - l1) quad[[3]] == l2 quad[[2]] + (1 - l2) quad[[4]],
        l1 >= 0 + convexThreshold,
        l1 <= 1 - convexThreshold,
        l2 >= 0 + convexThreshold,
        l2 <= 1 - convexThreshold
    }, {l1, l2}]] > 0;
intersectWithPlane[p1_, p2_, norm_, pointOnPlane_] := l1 p1 + (1 - l1) p2 /. Solve[ (l1 p1 + (1 - l1) p2).norm == pointOnPlane.norm, l1][[1]];

getEdgesFromFace[face_] := Table[face[[i+1]] - face[[Mod[i+1, Length[face]] + 1]], {i, 0, Length[face] - 1}];

quadri2DTests = {
    {
        {0, 0},
        {1, 0},
        {3/2, 7/2},
        {0, 3/2}
    },
    {
        {0, 0},
        {1, 0},
        {3, 3},
        {0, 1.5}
    },
    {
        {0, -1},
        {1.24453, 0},
        {0.325, 2.4532},
        {-0.234, 1.4335}
    },
    (* non convex *)
    {
        {0, 0},
        {-1, -3},
        {-1, 3},
        {3, 0}
    },

    {{0.47088352260634553, 0.38748337499380825}, {0.5191857247553109, 0.22942940224815433}, {0.8577913485838184, 0.08502211483622424}, {0.7504870010142595, 0.5340189184642878}},

    (* almost a triangle. Not "very" convex *)
    {{0.6699040278445817, 0.8436224721645679}, {0.8159470970115469, 0.12198382015814144}, {0.3212460475247576, 0.060176873941765394}, {0.05130863904584311, 0.0267112615210573}}
};

quadri2D = quadri2DTests[[2]];

convexCount = 0;
convexAttemptLimit = 1000;
For[attempt = 1, attempt <= convexAttemptLimit, attempt ++,
    quadri2D = Table[Random[], {i, 1, 4}, {j, 1, 2}];
    If[isConvex[quadri2D], convexCount += 1; Break[]];
];


(* If[attempt > convexAttemptLimit, 
Print["cannot construct convex random quadri"]; 
Print[{convexCount, attempt, convexCount/(attempt-1)//N, 25/36/3.0}];
Exit[]];
 *)

Print["quadrilateral is convex?"];
Print[isConvex[quadri2D]];
Print[quadri2D];

(* gridLimits = {-1/3, 1 + 1/3}; *)
gridLimits = {-1, 2};
gridGap = 1/3;
gridZ = -0.01;
grid3D = Table[
    {
        {{gridLimits[[1]], k, gridZ}, {gridLimits[[2]], k, gridZ}},
        {{k, gridLimits[[1]], gridZ}, {k, gridLimits[[2]], gridZ}}
    }, {k, gridLimits[[1]], gridLimits[[2]], gridGap}];
grid3D = Flatten[grid3D, 1];

quadri3D = Map[padZero, quadri2D];
innerCross = intersect[quadri2D[[1]], quadri2D[[3]], quadri2D[[2]], quadri2D[[4]]];
outerCross1 = intersect[quadri2D[[1]], quadri2D[[2]], quadri2D[[3]], quadri2D[[4]]];
outerCross2 = intersect[quadri2D[[1]], quadri2D[[4]], quadri2D[[2]], quadri2D[[3]]];

diagonal1 = intersect[outerCross1, outerCross2, quadri2D[[1]], quadri2D[[3]]];
diagonal2 = intersect[outerCross1, outerCross2, quadri2D[[2]], quadri2D[[4]]];

inner = padZero[innerCross];
o1 = padZero[outerCross1];
o2 = padZero[outerCross2];
d1 = padZero[diagonal1];
d2 = padZero[diagonal2];

c1 = (o1 + o2)/2;
c2 = (d1 + d2)/2;
r1 = Norm[(o2 - o1)/2];
r2 = Norm[(d2 - d1)/2];
solutionCenter = (c1 + c2)/2 + (r1^2 - r2^2)/Norm[c2 - c1]^2/2*(c2-c1);

solveResults = NSolve[{
        ({p1, p2, p3} - o1).({p1, p2, p3} - o2) == 0,
        ({p1, p2, p3} - d1).({p1, p2, p3} - d2) == 0,
        ({p1, p2, p3} - solutionCenter).(solutionCenter - inner) == 0,
        p3 >= 0
    }, {p1, p2, p3}, Reals];

If[Length[solveResults] == 0, 
    Print["No solution"]; 
    Exit[]
];

pSol = {p1, p2, p3}/.solveResults;

pSol = Sort[pSol, #1[[3]] > #2[[3]] &];

p = pSol[[1]];
Print[p];

planeNorm = Cross[o1-p, o2-p];
planeNorm = planeNorm/Norm[planeNorm];
Print["planeNorm"];
Print[planeNorm];

projectionMultiple = 1.5;
farthestVertex = MaximalBy[quadri3D, Abs[((#-p).planeNorm)]&][[1]];
projectedQ1 = (farthestVertex - p) * projectionMultiple + p;

projectedQuadri = Map[intersectWithPlane[p, #, planeNorm, projectedQ1]&, quadri3D];
projectedCenter = Mean[projectedQuadri];
planeNorm = planeNorm * If[(p - projectedCenter).planeNorm > 0, 1, -1];

orthPlaneNorm = RotationMatrix[Pi/2, {0, 0, 1}].planeNorm;
rotation3D = RotationMatrix[Pi/4, orthPlaneNorm];

expandFactor = 1.5;
expandedProjectedQuadri = Map[projectedCenter + (# - projectedCenter)*expandFactor &, projectedQuadri];

zs = Map[#[[3]]&, expandedProjectedQuadri];
zmin = Min[zs];
zmax = Max[zs];
otherTwoVertices = Select[expandedProjectedQuadri, #[[3]] > zmin + epsilon && #[[3]] < zmax - epsilon &];
xys = Map[#[[{1,2}]]&, otherTwoVertices];
screen = {
    {xys[[1]][[1]], xys[[1]][[2]], zmax},
    {xys[[2]][[1]], xys[[2]][[2]], zmax},
    {xys[[2]][[1]], xys[[2]][[2]], zmin},
    {xys[[1]][[1]], xys[[1]][[2]], zmin}
};

Print[Norm[projectedCenter - projectedQuadri[[1]]]];
If[Norm[projectedCenter - projectedQuadri[[1]]] > 1.5, 
    Print["Skip because the projected square is too big"];
    Exit[]
    ];

projectedEdges = getEdgesFromFace[projectedQuadri];

edgeLengths = Round[Join[
        Map[Norm, projectedEdges],
        Map[Norm, {projectedQuadri[[1]] - projectedQuadri[[3]], projectedQuadri[[2]] - projectedQuadri[[4]]}]/Sqrt[2]
        ], epsilon];

inners = Table[ Round[projectedEdges[[i + 1]].projectedEdges[[ Mod[ i + 1, 4] + 1  ]], epsilon] , {i, 0, 3}];
innerTally = Tally[inners];
edgeTally = Tally[edgeLengths];
(* Print[edgeTally]; *)
(* Print[innerTally]; *)
If[edgeTally[[1]][[2]] == 6 && innerTally[[1]][[2]] == 4, 
    Print["projectedquadri is a square"], 
    Print["projectedQuadri is not a square!!!!"]; Exit[]
];

(* rotate the whole thing by 45 deg to avoid visualization bug *)
p = rotation3D.p;
planeNorm = rotation3D.planeNorm;
projectedCenter = rotation3D.projectedCenter;
quadri3D = Map[(rotation3D.#) &, quadri3D];
grid3D = Map[(rotation3D.#) &, grid3D, {2}];
screen = Map[(rotation3D.#) &, screen];

range = 3;
lightingSource = {50, -50, 50};
shift = lightingSource - p;
lighting = {{"Point", White, lightingSource}};

thickness = 0.008;
startViewPoint = rotation3D.{0, 0, 1000};
endViewPoint = planeNorm * 1000;

lightDirection = (projectedCenter - p)/Norm[projectedCenter - p];
lightSideDirections = Orthogonalize[{lightDirection, {1, 0, 0}, {0, 0, 1}}][[{2, 3}]];

exportImage[viewAngle_, serialNumber_, frameFolder_] := Module[
    {
        i, v, dft, evRealImag, range, pointSize, graph, outputFolder, outputFileName, cropCommand
    },

    outputFolder = FileNameJoin[{imageFolder, frameFolder}];
    If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];

    outputFileName = "animation_frame_" <> IntegerString[serialNumber, 10, 4] <> ".png";
    outputFileName = FileNameJoin[{outputFolder, outputFileName}];

    Print["start and end viewpoint"];
    Print[startViewPoint];
    Print[endViewPoint];
    (* Print[startViewPoint.endViewPoint]; *)
    Print["viewpoint"];
    Print[startViewPoint Cos[viewAngle] + endViewPoint Sin[viewAngle]];
    graph1 = Graphics3D[
        {
            {
                FaceForm[Red], Polygon[Map[# + shift&, quadri3D]]
            },
            {
                FaceForm[Purple], Map[Cylinder[{#[[1]] + shift, #[[2]] + shift}, thickness]&, grid3D]
            },
            {
                FaceForm[Darker[Yellow, 0.5]], Cylinder[{p + shift - 0.1 lightDirection, p + shift - 0.2 lightDirection}, 0.05]
            },
            {
                FaceForm[Gray], Polygon[Map[# + shift&, screen]]
            }
        },
        ViewPoint -> startViewPoint Cos[viewAngle] + endViewPoint Sin[viewAngle],
        (* ViewVertical -> -planeNorm + {0, 0, 1}, *)
        Lighting -> lighting,
        Boxed -> False
    ];


    (* planeDirections = Orthogonalize[{planeNorm, {1, 0, 0}, {0, 0, 1}}][[{2, 3}]]; *)
(*     disk = ParametricPlot3D[r (Cos[t] planeDirections[[1]] + Sin[t] planeDirections[[2]]) + projectedCenter + shift, {r, 0, Norm[projectedQuadri[[1]] - projectedCenter] * expandFactor}, {t, 0, 2 Pi}, Boxed-> False, Mesh->None,
        PlotStyle -> Gray, PlotPoints -> 40
    ];
 *)

    cone = ParametricPlot3D[r (Cos[t] lightSideDirections[[1]] + Sin[t] lightSideDirections[[2]] + 1 * lightDirection ) - lightDirection * 0.15 + p + shift, {r, 0, 0.125}, {t, 0, 2 Pi}, Boxed-> False, Mesh->None,
        PlotStyle -> Yellow, PlotPoints -> 40
    ];

    graph = Show[graph1, cone];
(*     Print["plotrange"];
    originalPlotRange = FullOptions[graph, PlotRange];
    Print[originalPlotRange];
    newPlotRange = originalPlotRange * 0.9;
    Print[newPlotRange];
    graph = Show[graph1, disk, cone, PlotRange -> newPlotRange];
 *)
    cropCommand = "mogrify -crop " <> IntegerString[finalSize] <> "x" <> IntegerString[finalSize] <> "+" <> IntegerString[(imageSize[[1]]-finalSize)/2] <> "+" <> IntegerString[(imageSize[[2]]-finalSize)/2] <> " " <> outputFileName;

(* If you want the crop rectangle to start at top corner X: 50 Y: 100 and the crop rectangle to be of size W: 640 H:480, then use the command:

$ mogrify -crop 640x480+50+100 foo.png
 *)
    If[exportToPov,
        If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
        POVRayRender[graph, 
            Method -> "Triangle", 
            OutputPath -> outputFileName, 
            ImageSize -> imageSize
            , RenderOptions -> "+A0.001 -J +UA"
        ];
        Run[cropCommand];
        ,
        Export[outputFileName, graph, ImageSize -> finalSize];
    ];
    Print["Saved image to " <> outputFileName];
];

exportAnimation[] := Module[
    {
        shape, gifFileName, frameFolder, gifFolder, frameIndex, delay, gifCommand
    },

    shape = "quadri";
    gifFileName = shape 
        <> "_"
        <> IntegerString[frameCount]
        <> "_"
        <> IntegerString[RandomInteger[{10^10, 10^11-1}]];

    frameFolder = gifFileName <> "_frames";
    gifFolder = FileNameJoin[{imageFolder, shape}];
    outputFolder = FileNameJoin[{imageFolder, frameFolder}];
    If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];
    If[!DirectoryQ[gifFolder], CreateDirectory[gifFolder]];

    For[frameIndex = 1, frameIndex <= frameCount, frameIndex++,
        (* exportImage[Pi/2.0 * frameIndex/frameCount, frameIndex, frameFolder]; *)
        viewAngleWithoutHalfPi = (1 - Cos[frameIndex/frameCount * Pi 1.0])/2;
        Print["Frame #" <> IntegerString[frameIndex] <> " / " <> IntegerString[frameCount] <> " viewAngle=Pi/2* " <> ToString[viewAngleWithoutHalfPi]];
        exportImage[Pi/2.0 * viewAngleWithoutHalfPi, frameIndex, frameFolder];
    ];

    delay = If[frameCount >= 30, 3, Floor[500/frameCount]];
    gifCommand = "magick convert -duplicate 1,-2-1 -quiet +repage -alpha remove -loop 0 -delay " <> IntegerString[delay] <> " -dispose 2 " <> outputFolder <> "\\*animation*.png " <> gifFolder <> "\\" <> gifFileName <> ".gif";

    Print[gifCommand];
    Run[gifCommand];
    Print["Saved the animation to " <> gifFolder <> "\\" <> gifFileName <> ".gif"];

    mp4Command = "ffmpeg -i " <> gifFolder <> "\\" <> gifFileName <> ".gif" 
        <> " -movflags faststart -pix_fmt yuv420p -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" "
        <> gifFolder <> "\\" <> gifFileName <> ".mp4";

    mp4FileName = gifFolder <> "\\" <> gifFileName <> ".mp4";
    If[FileExistsQ[mp4FileName], DeleteFile[mp4FileName]];
    Print[mp4Command];
    Run[mp4Command];
    Print["Saved the video to " <> gifFolder <> "\\" <> gifFileName <> ".mp4"];
];

exportAnimation[];

(* exportImage[Pi/2 * (0/100), 0, "test"];
exportImage[Pi/2 * (1/100), 1, "test"];
exportImage[Pi/2 * (2/100), 2, "test"];
exportImage[Pi/2 * (50/100), 3, "test"];
exportImage[Pi/2 * (100/100), 4, "test"];
 *)