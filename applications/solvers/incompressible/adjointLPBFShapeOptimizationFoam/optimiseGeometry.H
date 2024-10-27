/* AM constraints algorithm
   author: jana.sasse.aachen@gmail.com */

Info << "Running geometric constraints algorithm" << endl;

/* -------------------------------------------------------------------- */
/* Use the linear density filter and design variables to get the
 intermediate densities (Equation 5.1) */

scalar sum_numerator = 0.0, sum_denominator = 0.0, rho_m;
label cellm;

forAll (rho_im, celli)
{
    forAll (allNeighbours[celli], m)
    {
        cellm = allNeighbours[celli][m];
        // alpha is scaled to [0, 1]
        rho_m = alpha[cellm]/alphaMax.value();        
        sum_numerator += weights[celli][m] * Vcell[cellm] * rho_m;
        sum_denominator += weights[celli][m] * Vcell[cellm];
    }

    rho_im[celli] = max(min(sum_numerator / sum_denominator, 1.0), 0.0);
    sum_numerator = 0.0, sum_denominator = 0.0;
}


/* -------------------------------------------------------------------- */
/* Use the nonlinear density filter and intermediate densities to get
 the physical densities (Equation 5.2) */

// eta: threshold parameter
scalar eta = 0.2;
// beta indicates how aggressive the filter is 
scalar beta = 0.2; 

forAll (rho_physical, celli)
{
    // NOTE [ZC20] instructs to include 0 and 1 as boundaries for rho_im,
    // but I found this leads to local errors, this one is more robust
    if (rho_im[celli] <= eta)
    {
        rho_physical[celli] = 
        (
            eta * 
            (
                Foam::exp(-beta * (1 - (rho_im[celli]/eta))) 
              - (1 - (rho_im[celli]/eta)) * Foam::exp(-beta)
            )
        );
    }
    else
    {
        rho_physical[celli] = 
        (
            (1-eta) * 
            (
                1 - Foam::exp(-beta * (rho_im[celli]-eta) / (1-eta)) 
              + ((rho_im[celli]-eta) * Foam::exp(-beta) / (1-eta))
            ) + eta
        );
    }
}


/* -------------------------------------------------------------------- */
/* Calculate geometric constraints using rho_physical */

// Overhang angle constraint (Equation 5.10)

scalarField xi_a(mesh.nCells(), 0.0);
scalarField xi_b(mesh.nCells(), 0.0);
scalarField xi_c(mesh.nCells(), 0.0);
scalarField cosAngle(mesh.nCells());
scalar cos_betaCrit = Foam::cos(0.78); // 45 degrees
label neigh;

forAll (alpha, celli)
{
    xCell = mesh.C()[celli].x();
    yCell = mesh.C()[celli].y();
    zCell = mesh.C()[celli].z();
    forAll (neighbour[celli], neighi)
    {
        neigh = neighbour[celli][neighi];
        xi_a[celli] += B_xRow[celli][neighi] * rho_physical[neigh];
        xi_b[celli] += B_yRow[celli][neighi] * rho_physical[neigh];
        xi_c[celli] += B_zRow[celli][neighi] * rho_physical[neigh];
    }
    // compute cosAngle (alpha in literature) = c/sqrt(a²+b²+c²)
    cosAngle[celli] = 
    (
        xi_a[celli] // we use xi_a here because we extrude in x-direction
      / (Foam::sqrt(sqr(xi_a[celli])+sqr(xi_b[celli])+sqr(xi_c[celli]))+SMALL)
    );

    // overhangCriterion < 0: these are the cells that are okay to be added
    // (based on the current structure) 
    overhangCriterion[celli] = 
    (
        (cosAngle[celli]-cos_betaCrit) 
      * Foam::sqrt(sqr(xi_a[celli]) + sqr(xi_b[celli]) + sqr(xi_c[celli]))
    );
}
    

// Horizontal minimum length and connectivity constraint (Equation 5.20)

scalarField rho_phase = rho_physical;
List<scalar> horLength(4, 0.0);
scalar horLength_1, horLength_2, horLength_3, horLength_4;
    
forAll (alpha, celli)
{
    horLength_1 = 0, horLength_2 = 0, horLength_3 = 0, horLength_4 = 0;

    // (Equation 5.19)
    if (rho_physical[celli] <= 0.5)
    {
        rho_phase[celli] = 1 - rho_physical[celli];
    }

    forAll (M_1[celli], i)
    {
        horLength_1 += (1 - rho_phase[M_1[celli][i]]);
    }
    // M_1[celli].size() is NOT always >= 1! (boundary region)
    horLength[0] = rho_phase[celli] - (horLength_1/max(1,M_1[celli].size()));
            
    forAll (M_2[celli], i)
    {
        horLength_2 += (1 - rho_phase[M_2[celli][i]]);
    }
    // M_2[celli].size() is NOT always >= 1! (boundary region)
    horLength[1] = rho_phase[celli] - (horLength_2/max(1,M_2[celli].size()));

    forAll (M_3[celli], i)
    {
        horLength_3 += (1 - rho_phase[M_3[celli][i]]);
    }
    // M_3[celli].size() is NOT always >= 1! (boundary region)
    horLength[2] = rho_phase[celli] - (horLength_3/max(1,M_3[celli].size()));

    forAll (M_4[celli], i)
    {
        horLength_4 += (1 - rho_phase[M_4[celli][i]]);
    }
    // M_4[celli].size() is NOT always >= 1! (boundary region)
    horLength[3] = rho_phase[celli] - (horLength_4/max(1,M_4[celli].size()));

    // minLengthCriterion should be smaller than 0 for at least one k! 
    minLengthCriterion[celli] = min(horLength);
}

/* -------------------------------------------------------------------- */
/* Update design variables */

optimizationProperties.lookup("lambdaU") >> lambdaU;
optimizationProperties.lookup("lambdaT") >> lambdaT;
optimizationProperties.lookup("lambdaC") >> lambdaC;
// compute individual sensitivities
volScalarField sensU = Ua & U;
volScalarField sensT = Ta * (T - average(T));
volScalarField sensC = Conca * (Conc - average(Conc));

// proposed design variable update
alphaDiff = 
(
    mesh.fieldRelaxationFactor("alpha")
  * (
        min
        (
            max
            (
                alpha + lambdaU*sensU + lambdaT*wT*sensT + lambdaC*wC*sensC
              , zeroAlpha
            )
          , alphaMax
        ) 
      - alpha
    )
);

// automatically reject all changes for scaffold structure
forAll(alpha, i)
{
    if (alphaFixed[i] == 1)
    {
        alpha[i] = alphaMax.value();
        alphaDiff[i] = 0;
    }
}

label proposedAddedCells = 0;
label proposedRemovedCells = 0;
label confirmedAddedCells = 0;
label confirmedRemovedCells = 0;
label sizeSolid = 0;

forAll(alpha, i)
{
    if (alphaDiff[i] > 0)
    {
        proposedAddedCells++;
    }
    if (alphaDiff[i] < 0)
    {
        proposedRemovedCells++;
    }
    if (alpha[i] > (0.5*alphaMax.value()))
    {
        sizeSolid++;
    }
}
Info<< "Number of solid cells before update: " << sizeSolid << endl;

zeroCells(alpha, inletCells);
//zeroCells(alpha, outletCells);

// Use geometric constraints to confirm or deny the proposed changes
forAll(alpha, i)
{
    if (alphaFixed[i] == 1)
    {
        alpha[i] = alphaMax.value();
        alphaDiff[i] = 0;
        continue;
    }

    // solid cell wants to stay solid, it is okay to build there
    // and it is happy there: confirm
    // fluid cell wants to become solid, it is okay to build there,
    // but it was also happy as fluid: confirm for now
    if 
    (
        (alphaDiff[i] > 0) 
     && (overhangCriterion[i] < 0) && (minLengthCriterion[i] <= 0)
    )
    {
        alpha[i] = min(alpha[i] + 5 * alphaDiff[i], alphaMax.value());
        confirmedAddedCells++;
        continue;
    }
    // fluid cell wants to become solid, it is okay to build there
    // and it was not happy as a fluid anyway: confirm
    if 
    (
        (alphaDiff[i] > 0) && (alpha[i] < 0.5*alphaMax.value())
     && (overhangCriterion[i] < 0) && (minLengthCriterion[i] > 0) 
    )
    {
        alpha[i] = min(alpha[i] + 5 * alphaDiff[i], alphaMax.value());
        confirmedAddedCells++;
        continue;
    }
    // solid cell wants to become fluid, although AM restrictions
    // would be fine with it staying solid: confirm for now
    // fluid cell wants to stay fluid and AM restrictions don't care: confirm
    if 
    (
        (alphaDiff[i] < 0) 
     && (overhangCriterion[i] < 0) && (minLengthCriterion[i] <= 0)
    )
    {
        alpha[i] = max(alpha[i] + 5 * alphaDiff[i], zeroAlpha.value());
        confirmedRemovedCells++;
        continue;
    }
    // solid cell wants to become fluid, overhang was fine there
    // but it was not happy in hor. length: confirm
    if 
    (
        (alphaDiff[i] < 0) && (alpha[i] > 0.5*alphaMax.value())
     && (overhangCriterion[i] < 0) && (minLengthCriterion[i] > 0)
    )
    {
        alpha[i] = max(alpha[i] + 5 * alphaDiff[i], zeroAlpha.value());
        confirmedRemovedCells++;
        continue;
    }
    // solid cell wants to become fluid, overhang was not fine anyway:
    // confirm for now
    // fluid cell wants to stay fluid and AM restrictions are happy
    // this way: confirm
    if 
    (
        (alphaDiff[i] < 0) 
     && (overhangCriterion[i] >= 0) && (minLengthCriterion[i] <= 0)
    )
    {
        alpha[i] = max(alpha[i] + 5 * alphaDiff[i], zeroAlpha.value());
        confirmedRemovedCells++;
        continue;
    }
    // solid cell wants to become fluid, overhang was not fine there
    // and it was not happy in hor. length: confirm
    if 
    (
        (alphaDiff[i] < 0) && (alpha[i] > 0.5*alphaMax.value())
     && (overhangCriterion[i] >= 0) && (minLengthCriterion[i] > 0)
    )
    {
        alpha[i] = max(alpha[i] + 5 * alphaDiff[i], zeroAlpha.value());
        confirmedRemovedCells++;
        continue;
    }
    // reject all other changes automatically and continue
}

// print final number of solid cells after application of constraints
sizeSolid = 0;
forAll(alpha, i)
{
    if (alpha[i] > (0.5*alphaMax.value()))
    {
        sizeSolid++;
    }
}
label rejectedAddedCells = proposedAddedCells - confirmedAddedCells;
label rejectedRemovedCells = proposedRemovedCells - confirmedRemovedCells;

Info<< "    " 
    << proposedAddedCells << " cells want to become solid. It was confirmed for " 
    << confirmedAddedCells << " cells and rejected for " 
    << rejectedAddedCells << " cells." << endl;
Info<< "    " 
    << proposedRemovedCells << " cells want to become fluid. It was confirmed for "
    << confirmedRemovedCells << " cells and rejected for " 
    << rejectedRemovedCells << " cells." << endl;
Info << "Number of solid cells after update: " << sizeSolid << "\n" << endl;
