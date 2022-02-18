(* ::Package:: *)

BeginPackage[ "PQSHelper`" ];

AnalyticSolve::usage = 
	"AnalyticSolve[eq, n, x, \[Sigma], c] computes the analytic estimated solution to eq == 0.";

Begin[ "`Private`" ]

AnalyticSolve[eq1_, n_, x_, \[Sigma]_, c_] := 
	Module[ {eq2, eq3, list1, list2, list3, notrig, htrigident, trigident, order},
		notrig = {
			Cosh[_] -> cosh, Sinh[_]->sinh,
			Tanh[_] -> sinh / cosh, Coth[_] -> cosh / sinh, 
			Sech[_] -> 1/cosh,  Csch[_] -> 1/sinh
		};
		htrigident = {
			cosh^k_-> (1+sinh^2)^(Floor[k/2]) cosh^(Mod[k,2])
		};
		trigident = {
			Cos[x \[Sigma]]^k_-> (1-Sin[x \[Sigma]]^2)^(Floor[k/2]) Cos[x \[Sigma]]^(Mod[k,2])
		};
		eq2 = eq1 //. notrig // Together;
		eq3 = ((eq2 // Numerator) //. htrigident) // Expand;
		list1 = CoefficientList[eq3, {sinh, cosh}];
		(* For u1, we had denom = 2*7 = 14, so numerator orders 11 *)
		(* For u2, we had denom = 2*9 = 18, so numerator orders 13 *)
		(* For u3, we had denom = 2*15 = 30, so numerator orders 23 *)
		order = ((Length[CoefficientList[#, cosh]] - 1) + (Length[CoefficientList[#, sinh]] - 1))&[eq2 // Denominator];
		list2 = {
			list1[[order - (2 n + 1), 2]],
			list1[[order - (2 n + 1) + 1, 1]]
		};
		list3 = CoefficientList[list2 /. trigident, {Cos[x \[Sigma]], Sin[x \[Sigma]]}] // Flatten // Select[Not@*PossibleZeroQ];
		Solve[list3 == 0, (i |-> c[i]) /@ Range[0, 2 * n + 1]][[1]] // FullSimplify
	]
	
End[];

EndPackage[];
