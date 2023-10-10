/* PetscErrorCode GetDikeContr(JacRes *jr,                                                                                                                                
                            PetscScalar *phRat,          // phase ratios in the control volume   
                            PetscInt &AirPhase,                                                                           
                            PetscScalar &dikeRHS,
                            PetscScalar &y_c,
                            PetscInt J,
                            PetscInt I,
                            PetscScalar sxx_eff_ave_cell) 
                                             
{
  BCCtx       *bc;
  Dike        *dike;
  Ph_trans_t  *CurrPhTr;
  PetscScalar v_spread, M, left, right, front, back;
  PetscScalar y_distance, tempdikeRHS;
  PetscInt    i, nD, nPtr, numDike, numPhtr, nsegs;

  // *djking
  PetscScalar P_comp, P_star, Pm;
  PetscScalar div_max, M_val, dike_or, zeta;
  PetscInt    L, sy, sx;
  PetscMPIInt rank;
  //

  PetscFunctionBeginUser;

// *revisit (can this stuff be called from jr like in Dike_k_heatsource?? if not use dikeContr_temp to fix)
  numDike    = jr->dbdike->numDike; // number of dikes
  numPhtr    = jr->dbm->numPhtr;
  bc         = jr->bc;

  nPtr = 0;
  nD = 0;

    for(nPtr=0; nPtr<numPhtr; nPtr++)   // loop over all phase transitions blocks
    {
      // access the parameters of the phasetranstion block
      CurrPhTr = jr->dbm->matPhtr+nPtr;
      
      for(nD = 0; nD < numDike; nD++) // loop through all dike blocks
      {
          // access the parameters of the dike depending on the dike block
          dike = jr->dbdike->matDike+nD;
	  
	        // access the phase ID of the dike parameters of each dike
          i = dike->PhaseID;
	  
          if(CurrPhTr->ID == dike->PhaseTransID)  // compare the phaseTransID associated with the dike with the actual ID of the phase transition in this cell           
          {
//      PetscPrintf(PETSC_COMM_WORLD,"if CurrPhTr \n");
//   if(phRat[i]>0)
//    {
//      PetscPrintf(PETSC_COMM_WORLD,"phRat[i] = %g, ", phRat[i]);
//    }
	           // check if the phase ratio of a dike phase is greater than 0 in the current cell
	           if(phRat[i]>0 && CurrPhTr->celly_xboundR[J] > CurrPhTr->celly_xboundL[J])
		         {
                nsegs=CurrPhTr->nsegs;
//        PetscPrintf(PETSC_COMM_WORLD,"phRat > 0 && xboundR > xboundL, ");

                P_comp = - sxx_eff_ave_cell * 1e9 + dike->Ts; // *revisit (scale; multiplied by 1e9 to convert to Pa)

//        PetscPrintf(PETSC_COMM_WORLD,"sxx_eff_ave_cell = %g, ", sxx_eff_ave_cell);
//        PetscPrintf(PETSC_COMM_WORLD,"P_comp = %g, ", P_comp);

		            if(dike->Mb == dike->Mf && dike->Mc < 0.0) // spatially constant M
		            {
//        PetscPrintf(PETSC_COMM_WORLD,"spatially constant M: ");
		              M = dike->Mf;
		              v_spread = PetscAbs(bc->velin);
		              left = CurrPhTr->celly_xboundL[J];
		              right = CurrPhTr->celly_xboundR[J];
		              
                  if(jr->ctrl.var_M) // *djking
                  {
//        PetscPrintf(PETSC_COMM_WORLD,"var_M on\n");
                    M_val = M * PetscSqrtReal(PetscAbs(-P_comp / (dike->knee + PetscAbs(P_comp))));
//                    div_max = M_val * 2 * v_spread / PetscAbs(left-right); // maximum divergence allowed
//                    dike_or = M * 2 * v_spread * 10; // dike opening rate (km/Myr) *revisit (scale??)
                    div_max = M_val * 2 * (v_spread * 315.57599999999996 / 100 / 365.25 / 24 / 60/ 60)/ (PetscAbs(left-right) * 1000); //*hardcoded scale
                    dike_or = M * 2 * (v_spread * 315.57599999999996 / 100 / 365.25 / 24 / 60/ 60); //*hardcoded scale
                    
        PetscPrintf(PETSC_COMM_WORLD,"P_comp = %g (MPa), ", P_comp/1e6);
        PetscPrintf(PETSC_COMM_WORLD,"sxx_eff_ave_cell = %g (MPa), ", sxx_eff_ave_cell*1e3);
        PetscPrintf(PETSC_COMM_WORLD,"div_max = %g, ", div_max);
        PetscPrintf(PETSC_COMM_WORLD,"M_val = %g, ", M_val);
        PetscPrintf(PETSC_COMM_WORLD,"v_spread = %g, ", v_spread * 315.57599999999996 / 100 / 365.25 / 24 / 60/ 60); //*hardcoded scale
//        PetscPrintf(PETSC_COMM_WORLD,"v_spread = %g, ", v_spread);
        PetscPrintf(PETSC_COMM_WORLD,"left = %g, ", left);
        PetscPrintf(PETSC_COMM_WORLD,"right = %g \n", right);

                    if(P_comp < 0) // diking occurs
                    {
                      zeta = -(dike->A * dike->zeta_0 / (P_comp + dike->B) + P_comp / div_max);
                      tempdikeRHS = - P_comp / zeta; // *revisit (scale;  non-dim time and lengths...)

        PetscPrintf(PETSC_COMM_WORLD,"diking --> ");
//         PetscPrintf(PETSC_COMM_WORLD,"M = %g, ", M_val);
        PetscPrintf(PETSC_COMM_WORLD,"tempdikeRHS = %g \n", tempdikeRHS);
//         PetscPrintf(PETSC_COMM_WORLD,"div_max = %g, ", div_max);
//        PetscPrintf(PETSC_COMM_WORLD,"left = %g, ", left);
//        PetscPrintf(PETSC_COMM_WORLD,"right = %g \n", div_max);
//        PetscPrintf(PETSC_COMM_WORLD,"phRat[i] = %g, ", phRat[i]);
//        PetscPrintf(PETSC_COMM_WORLD,"sxx_eff_ave_cell = %g, ", sxx_eff_ave_cell);
//        PetscPrintf(PETSC_COMM_WORLD,"div_max = %g, ", div_max);
        PetscPrintf(PETSC_COMM_WORLD,"P_comp = %g, ", P_comp);
        PetscPrintf(PETSC_COMM_WORLD,"zeta = %g \n", zeta);
                    }
                    else // diking DOES NOT occur
                    {
                      tempdikeRHS = 0.0;
//        PetscPrintf(PETSC_COMM_WORLD,"not diking --> ");
//        PetscPrintf(PETSC_COMM_WORLD,"M = %g \n", tempdikeRHS);
                    }

//        PetscPrintf(PETSC_COMM_WORLD,"M = %g \n", tempdikeRHS);

                  }
                  else
                  {
                    tempdikeRHS = M * 2 * v_spread / PetscAbs(left-right);
        PetscPrintf(PETSC_COMM_WORLD,"var_M off\n");
        PetscPrintf(PETSC_COMM_WORLD,"M = %g, ", M);
        PetscPrintf(PETSC_COMM_WORLD,"tempdikeRHS = %g \n", tempdikeRHS);
                  }


		            }

		            else if(dike->Mc >= 0.0)   // Mf, Mc and Mb
		            {
                  left = CurrPhTr->celly_xboundL[J];
                  right = CurrPhTr->celly_xboundR[J];
                  front = CurrPhTr->ybounds[0];
                  back = CurrPhTr->ybounds[2*nsegs-1];
                  v_spread = PetscAbs(bc->velin);

                  if(y_c >= dike->y_Mc)
                  {
                      // linear interpolation between different M values, Mc is M in the middle, acts as M in front, Mb is M in back 
                      y_distance = y_c - dike->y_Mc;
                      M = dike->Mc + (dike->Mb - dike->Mc) * (y_distance / (back - dike->y_Mc));
                      //tempdikeRHS = M * 2 * v_spread / PetscAbs(left - right);
                  }
                  
                  else
                  {
                      // linear interpolation between different M values, Mf is M in front, Mc acts as M in back  
                      y_distance = y_c - front;
                      M = dike->Mf + (dike->Mc - dike->Mf) * (y_distance / (dike->y_Mc - front));
                      //tempdikeRHS = M * 2 * v_spread / PetscAbs(left - right);
                  }
		            }
		            else if(dike->Mb != dike->Mf && dike->Mc < 0.0)   // only Mf and Mb, they are different
		            {
                  left = CurrPhTr->celly_xboundL[J];
                  right = CurrPhTr->celly_xboundR[J];
                  back = CurrPhTr->ybounds[2*nsegs-1];

                  v_spread = PetscAbs(bc->velin);
        
                  // linear interpolation between different M values, Mf is M in front, Mb is M in back
                  y_distance = y_c - front;
                  M = dike->Mf + (dike->Mb - dike->Mf) * (y_distance / (back - front));
                  //tempdikeRHS = M * 2 * v_spread / PetscAbs(left - right);
		            }
		            else   // Mb and Mf don't exist (which should not occurr)
		            {
		              SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "No values for Mb and Mf. Dike option invalid!");
		            }

                // Divergence
		              dikeRHS += (phRat[i]+phRat[AirPhase])*tempdikeRHS;  // Give full divergence if cell is part dike part air

		        }  //close if phRat and xboundR>xboundL  
	        }  // close phase transition and dike phase ID comparison 
	    }  // close dike block loop
    }  // close phase transition block loop
  
  PetscFunctionReturn(0);
} */