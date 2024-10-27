/* AM constraints algorithm
   author: jana.sasse.aachen@gmail.com */

Info << "Precomputations for geometric constraints algorithm" << endl;

/* -------------------------------------------------------------------- */
scalar yNeigh, zNeigh, xCell, yCell, zCell;
scalar yNeighNeigh, zNeighNeigh, radius, distance, init;
label n1, n2, n3, n4, n5, n6;

// list for all neighbours of all cells
DynamicList<List<label>> neighbour(mesh.nCells());
DynamicList<List<label>> allNeighbours(mesh.nCells());
DynamicList<List<scalar>> B_xRow(mesh.nCells());
DynamicList<List<scalar>> B_yRow(mesh.nCells());
DynamicList<List<scalar>> B_zRow(mesh.nCells());
DynamicList<List<label>> M_1(mesh.nCells());
DynamicList<List<label>> M_2(mesh.nCells());
DynamicList<List<label>> M_3(mesh.nCells());
DynamicList<List<label>> M_4(mesh.nCells());
DynamicList<List<label>> sameLayer(mesh.nCells());

const scalarField&      Vcell  = mesh.V(); 
DynamicList<List<scalar>> weights(mesh.nCells());

forAll (alpha, celli)
{
    /* -------- determine neighbours and weights of celli -------- */
    DynamicList<label> N_i_sameLayer;
    DynamicList<label> N_i_lowerLayer;
    DynamicList<label> N_i_upperLayer;
    DynamicList<scalar> w_i_sameLayer;
    DynamicList<scalar> w_i_lowerLayer;
    DynamicList<scalar> w_i_upperLayer;
    xCell = mesh.C()[celli].x();
    yCell = mesh.C()[celli].y();
    zCell = mesh.C()[celli].z();

    // filter radius "R"
    radius = 0.2 * Foam::sqrt(Foam::sqrt(Vcell[celli]));

    // find points corresponding to celli
    const labelList& cellPoints = mesh.cellPoints()[celli];
    forAll (cellPoints, pointi)
    {
        // identify all cells that use this point
        const labelList& pointCells = mesh.pointCells()[cellPoints[pointi]];
        forAll (pointCells, pointCelli)
        {
            // check for each cell using this point if ...
            // (a) ... it is not already in our list of neighbours of celli
            // (b) ... this cell is not celli itself
            // (c) ... cell is in same layer as celli
            if (
                (findIndex(N_i_sameLayer, pointCells[pointCelli]) == -1)
             && (pointCells[pointCelli] != celli)
             && (mesh.C()[pointCells[pointCelli]].x() > (xCell - SMALL))
             && (mesh.C()[pointCells[pointCelli]].x() < (xCell + SMALL))
            )
            {
                N_i_sameLayer.append(pointCells[pointCelli]);
                distance = Foam::sqrt
                    (
                        sqr(mesh.C()[pointCells[pointCelli]].x()-xCell) 
                    + sqr(mesh.C()[pointCells[pointCelli]].y()-yCell)
                    + sqr(mesh.C()[pointCells[pointCelli]].z()-zCell)
                    );
                w_i_sameLayer.append(radius-distance);
            }

            // check for each cell using this point if ...
            // (a) ... it is not already in our list of neighbours of celli
            // (b) ... this cell is not celli itself
            // (c) ... cell is in the layer below celli
            if (
                (findIndex(N_i_lowerLayer, pointCells[pointCelli]) == -1)
             && (pointCells[pointCelli] != celli)
             && (mesh.C()[pointCells[pointCelli]].x() < (xCell - SMALL))
            )
            {
                N_i_lowerLayer.append(pointCells[pointCelli]);
                distance = Foam::sqrt
                    (
                        sqr(mesh.C()[pointCells[pointCelli]].x()-xCell) 
                    + sqr(mesh.C()[pointCells[pointCelli]].y()-yCell)
                    + sqr(mesh.C()[pointCells[pointCelli]].z()-zCell)
                    );
                w_i_lowerLayer.append(radius-distance);
            }

            // check for each cell using this point if ...
            // (a) ... it is not already in our list of neighbours of celli
            // (b) ... this cell is not celli itself
            // (c) ... cell is in the layer above celli
            if (
                (findIndex(N_i_upperLayer, pointCells[pointCelli]) == -1)
             && (pointCells[pointCelli] != celli)
             && (mesh.C()[pointCells[pointCelli]].x() > (xCell + SMALL))
            )
            {
                N_i_upperLayer.append(pointCells[pointCelli]);
                distance = Foam::sqrt
                    (
                        sqr(mesh.C()[pointCells[pointCelli]].x()-xCell) 
                    + sqr(mesh.C()[pointCells[pointCelli]].y()-yCell)
                    + sqr(mesh.C()[pointCells[pointCelli]].z()-zCell)
                    );
                w_i_upperLayer.append(radius-distance);
            }
        }

    }

    // NOTE: we now have to remain in the order defined by sameLayer and lowerLayer
    // first cell in neighbour is i
    neighbour[celli].append(celli);
    neighbour[celli].append(N_i_sameLayer);
    neighbour[celli].append(N_i_lowerLayer);
    allNeighbours[celli].append(celli);
    allNeighbours[celli].append(N_i_sameLayer);
    allNeighbours[celli].append(N_i_lowerLayer);
    allNeighbours[celli].append(N_i_upperLayer);
    weights[celli].append(0.0);
    weights[celli].append(w_i_sameLayer);
    weights[celli].append(w_i_lowerLayer);
    weights[celli].append(w_i_upperLayer);

    // need this for Part II of minLengthConstraint precomputations
    sameLayer[celli].append(N_i_sameLayer);

    /* ---- Assembly of B_i for (in-plane) non-cartesian mesh ---- */

    // assign B_xRow
    // n1 = number of cells that get + Label (add 1 to include celli)
    n1 = N_i_sameLayer.size() + 1; 
    init = 1.0/n1;
    List<scalar> B_xRow_sameLayer(n1, init);
    B_xRow[celli].append(B_xRow_sameLayer);
    // n2 = number of cells that get - Label
    n2 = N_i_lowerLayer.size();
    if (n2 > 0)
    {
        init = -1.0/n2;
        List<scalar> B_xRow_lowerLayer(n2, init);
        B_xRow[celli].append(B_xRow_lowerLayer);
    }

    // assign B_yRow, B_zRow
    n3=0, n4=0, n5=0, n6=0;
    init = 0.0;
    List<scalar> B_yRow_sameLayer(n1-1, init);
    List<scalar> B_yRow_lowerLayer(n2, init);
    List<scalar> B_zRow_sameLayer(n1-1, init);
    List<scalar> B_zRow_lowerLayer(n2, init);
    forAll (N_i_sameLayer, neighi)
    {
        yNeigh = mesh.C()[N_i_sameLayer[neighi]].y();
        zNeigh = mesh.C()[N_i_sameLayer[neighi]].z();
        
        if (yNeigh > (yCell + SMALL))
        {
            n3++;
            B_yRow_sameLayer[neighi] = 1.0;
        }
        if (yNeigh < (yCell - SMALL))
        {
            n4++;
            B_yRow_sameLayer[neighi] = -1.0;
        }
        if (zNeigh > (zCell + SMALL))
        {
            n5++;
            B_zRow_sameLayer[neighi] = 1.0;
        }
        if (zNeigh < (zCell - SMALL))
        {
            n6++;
            B_zRow_sameLayer[neighi] = -1.0;
        }

        //precomputations for minLengthConstraint PART I
        const labelList& faceNeighs = mesh.cellCells()[celli];
        const labelList& neighFaceNeighs 
            = mesh.cellCells()[N_i_sameLayer[neighi]];
        if ((yNeigh > (yCell - SMALL)) && (zNeigh < (zCell + SMALL)))
        {
            if(findIndex(M_1[celli], N_i_sameLayer[neighi]) == -1)
            {
                M_1[celli].append(N_i_sameLayer[neighi]);
            }
            forAll(faceNeighs, faceNeighi)
            {
                // if the cell is in faceNeighs list of both celli and neighi
                // AND the cell has not been added to M yet
                if
                (
                 (findIndex(neighFaceNeighs, faceNeighs[faceNeighi]) != -1)
                 && (findIndex(M_1[celli], faceNeighs[faceNeighi]) == -1)
                )
                {
                    M_1[celli].append(faceNeighs[faceNeighi]);
                }
            }
        }
        if ((yNeigh < (yCell + SMALL)) && (zNeigh < (zCell + SMALL)))
        {
            if(findIndex(M_2[celli], N_i_sameLayer[neighi]) == -1)
            {
                M_2[celli].append(N_i_sameLayer[neighi]);
            }
            forAll(faceNeighs, faceNeighi)
            {
                // if the cell is in faceNeighs list of both celli and neighi
                // AND the cell has not been added to M yet
                if
                (
                 (findIndex(neighFaceNeighs, faceNeighs[faceNeighi]) != -1)
                 && (findIndex(M_2[celli], faceNeighs[faceNeighi]) == -1)
                )
                {
                    M_2[celli].append(faceNeighs[faceNeighi]);
                }
            }
        }
        if ((yNeigh < (yCell + SMALL)) && (zNeigh > (zCell - SMALL)))
        {
            if(findIndex(M_3[celli], N_i_sameLayer[neighi]) == -1)
            {
                M_3[celli].append(N_i_sameLayer[neighi]);
            }
            forAll(faceNeighs, faceNeighi)
            {
                // if the cell is in faceNeighs list of both celli and neighi
                // AND the cell has not been added to M yet
                if
                (
                 (findIndex(neighFaceNeighs, faceNeighs[faceNeighi]) != -1)
                 && (findIndex(M_3[celli], faceNeighs[faceNeighi]) == -1)
                )
                {
                    M_3[celli].append(faceNeighs[faceNeighi]);
                }
            }
        }
        if ((yNeigh > (yCell - SMALL)) && (zNeigh > (zCell - SMALL)))
        {
            if(findIndex(M_4[celli], N_i_sameLayer[neighi]) == -1)
            {
                M_4[celli].append(N_i_sameLayer[neighi]);
            }
            forAll(faceNeighs, faceNeighi)
            {
                // if cell is in faceNeighs list of both celli and neighi
                // AND the cell has not been added to M yet
                if
                (
                 (findIndex(neighFaceNeighs, faceNeighs[faceNeighi]) != -1)
                 && (findIndex(M_4[celli], faceNeighs[faceNeighi]) == -1)
                )
                {
                    M_4[celli].append(faceNeighs[faceNeighi]);
                }
            }
        }
    }

    forAll (N_i_lowerLayer, neighi)
    {
        yNeigh = mesh.C()[N_i_lowerLayer[neighi]].y();
        zNeigh = mesh.C()[N_i_lowerLayer[neighi]].z();

        if (yNeigh > (yCell + SMALL))
        {
            n3++;
            B_yRow_lowerLayer[neighi] = 1.0;
        }
        if (yNeigh < (yCell - SMALL))
        {
            n4++;
            B_yRow_lowerLayer[neighi] = -1.0;
        }
        if (zNeigh > (zCell + SMALL))
        {
            n5++;
            B_zRow_lowerLayer[neighi] = 1.0;
        }
        if (zNeigh < (zCell - SMALL))
        {
            n6++;
            B_zRow_lowerLayer[neighi] = -1.0;
        }
    }

    // every cell found for n3 gets assigned +0.5/n3 for B_yRow, 
    // and every cell found for n4 gets assigned -0.5/n4 for B_yRow 
    // every cell found for n5 gets assigned +0.5/n5 for B_zRow, 
    // and every cell found for n6 gets assigned -0.5/n6 for B_zRow
    B_yRow[celli].append(B_yRow_sameLayer);
    B_yRow[celli].append(B_yRow_lowerLayer);
    forAll (B_yRow[celli], neighi)
    {
        if (B_yRow[celli][neighi] > 0)
        {
            B_yRow[celli][neighi] *= 0.5/n3;
        }
        if (B_yRow[celli][neighi] < 0)
        {
            B_yRow[celli][neighi] *= 0.5/n4;
        }
    }
    B_zRow[celli].append(B_zRow_sameLayer);
    B_zRow[celli].append(B_zRow_lowerLayer);
    forAll (B_zRow[celli], neighi)
    {
        if (B_zRow[celli][neighi] > 0)
        {
            B_zRow[celli][neighi] *= 0.5/n5;
        }
        if (B_zRow[celli][neighi] < 0)
        {
            B_zRow[celli][neighi] *= 0.5/n6;
        }
    }

}

// precomputations for minLengthConstraint PART II;
// already found first-level neighs in loop above,
// now we need to append with second-level neighbours here
// (can't do it earlier because this requires all neighbours)
label currentNeigh, currentNeighNeigh;
forAll (alpha, celli)
{

    yCell = mesh.C()[celli].y();
    zCell = mesh.C()[celli].z();
    DynamicList<label> tmp1, tmp2, tmp3, tmp4;
    forAll (M_1[celli], neighi)
    {
        currentNeigh = M_1[celli][neighi];
        yNeigh = mesh.C()[currentNeigh].y();
        zNeigh = mesh.C()[currentNeigh].z();
        // these are first degree neighbours, now find their neighbours
        // and put them in some tmp-structure to append to M_1 later
        forAll (sameLayer[currentNeigh], neighNeighi)
        {
            currentNeighNeigh = sameLayer[currentNeigh][neighNeighi];
            yNeighNeigh = mesh.C()[currentNeighNeigh].y();
            zNeighNeigh = mesh.C()[currentNeighNeigh].z();
            // criteria 2, 3 and 4 remove celli and any duplicates
            if (((yNeighNeigh > (yNeigh-SMALL)) && (zNeighNeigh < (zNeigh+SMALL)))
                && (currentNeighNeigh != celli)
                && (findIndex(tmp1, currentNeighNeigh) == -1)
                && (findIndex(M_1[celli], currentNeighNeigh) == -1))
            {
                tmp1.append(currentNeighNeigh);
            }
        }   
    }
    M_1[celli].append(tmp1);

    forAll (M_2[celli], neighi)
    {
        currentNeigh = M_2[celli][neighi];
        yNeigh = mesh.C()[currentNeigh].y();
        zNeigh = mesh.C()[currentNeigh].z();
        // these are first degree neighbours, now find their neighbours
        // and put them in some tmp-structure to append to M_2 later
        forAll (sameLayer[currentNeigh], neighNeighi)
        {
            currentNeighNeigh = sameLayer[currentNeigh][neighNeighi];
            yNeighNeigh = mesh.C()[currentNeighNeigh].y();
            zNeighNeigh = mesh.C()[currentNeighNeigh].z();
            // criteria 2, 3 and 4 remove celli and any duplicates
            if (((yNeighNeigh < (yNeigh+SMALL)) && (zNeighNeigh < (zNeigh+SMALL)))
                && (currentNeighNeigh != celli)
                && (findIndex(tmp2, currentNeighNeigh) == -1)
                && (findIndex(M_2[celli], currentNeighNeigh) == -1))
            {
                tmp2.append(currentNeighNeigh);
            }
        }   
    }
    M_2[celli].append(tmp2);

    forAll (M_3[celli], neighi)
    {
        currentNeigh = M_3[celli][neighi];
        yNeigh = mesh.C()[currentNeigh].y();
        zNeigh = mesh.C()[currentNeigh].z();
        // these are first degree neighbours, now find their neighbours
        // and put them in some tmp-structure to append to M_3 later
        forAll (sameLayer[currentNeigh], neighNeighi)
        {
            currentNeighNeigh = sameLayer[currentNeigh][neighNeighi];
            yNeighNeigh = mesh.C()[currentNeighNeigh].y();
            zNeighNeigh = mesh.C()[currentNeighNeigh].z();
            // criteria 2, 3 and 4 remove celli and any duplicates
            if (((yNeighNeigh < (yNeigh+SMALL)) && (zNeighNeigh > (zNeigh-SMALL)))
                && (currentNeighNeigh != celli)
                && (findIndex(tmp3, currentNeighNeigh) == -1)
                && (findIndex(M_3[celli], currentNeighNeigh) == -1))
            {
                tmp3.append(currentNeighNeigh);
            }
        }   
    }
    M_3[celli].append(tmp3);

    forAll (M_4[celli], neighi)
    {
        currentNeigh = M_4[celli][neighi];
        yNeigh = mesh.C()[currentNeigh].y();
        zNeigh = mesh.C()[currentNeigh].z();
        // these are first degree neighbours, now find their neighbours
        // and put them in some tmp-structure to append to M_4 later
        forAll (sameLayer[currentNeigh], neighNeighi)
        {
            currentNeighNeigh = sameLayer[currentNeigh][neighNeighi];
            yNeighNeigh = mesh.C()[currentNeighNeigh].y();
            zNeighNeigh = mesh.C()[currentNeighNeigh].z();
            // criteria 2, 3 and 4 remove celli and any duplicates
            if (((yNeighNeigh > (yNeigh-SMALL)) && (zNeighNeigh > (zNeigh-SMALL)))
                && (currentNeighNeigh != celli)
                && (findIndex(tmp4, currentNeighNeigh) == -1)
                && (findIndex(M_4[celli], currentNeighNeigh) == -1))
            {
                tmp4.append(currentNeighNeigh);
            }
        }   
    }
    M_4[celli].append(tmp4);
    
}
