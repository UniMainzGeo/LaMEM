/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This software was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   solVar.h
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//.....................   FORWARD TYPE DECLARATIONS   .......................
//---------------------------------------------------------------------------
#ifndef __types_h__
#define __types_h__
//---------------------------------------------------------------------------

typedef struct Scaling        Scaling       ;
typedef struct TSSol          TSSol         ;
typedef struct DBMat          DBMat         ;
typedef struct FDSTAG         FDSTAG        ;
typedef struct AdvCtx         AdvCtx        ;
typedef struct BCCtx          BCCtx         ;
typedef struct JacRes         JacRes        ;
typedef struct PVAVD          PVAVD         ;
typedef struct PVOut          PVOut         ;
typedef struct PVMark         PVMark        ;
typedef struct PVSurf         PVSurf        ;
typedef struct FB             FB            ;
typedef struct FreeSurf       FreeSurf      ;

typedef struct Marker         Marker        ;
typedef struct SolVarDev      SolVarDev     ;
typedef struct SolVarBulk     SolVarBulk    ;
typedef struct SolVarCell     SolVarCell    ;
typedef struct SolVarEdge     SolVarEdge    ;
typedef struct Soft_t         Soft_t        ;
typedef struct Material_t     Material_t    ;
typedef struct MatParLim      MatParLim     ;
typedef struct Tensor2RS      Tensor2RS     ;
typedef struct Tensor2RN      Tensor2RN     ;
typedef struct Pair           Pair          ;


typedef struct NLSol          NLSol         ;
typedef struct VelInterp      VelInterp     ;
typedef struct AdvVelCtx      AdvVelCtx     ;
typedef struct InterpFlags    InterpFlags   ;
typedef struct DOFIndex       DOFIndex      ;

//typedef struct ConstEqCtx     ConstEqCtx    ;
//typedef struct Polygon2D      Polygon2D     ;
//typedef struct ModParam       ModParam      ;
//typedef struct ObjFunct       ObjFunct      ;
//typedef struct GravitySurvey  GravitySurvey ;
//typedef struct MeshSeg1D      MeshSeg1D     ;
//typedef struct Discret1D      Discret1D     ;
//typedef struct NumCorner      NumCorner     ;
//typedef struct BCBlock        BCBlock       ;
//typedef struct DBox           DBox          ;
//typedef struct OutBuf         OutBuf        ;
//typedef struct OutVec         OutVec        ;
//typedef struct OutMask        OutMask       ;
//typedef struct GeomPrim       GeomPrim      ;
//typedef struct AVDCell        AVDCell       ;
//typedef struct AVDChain       AVDChain      ;
//typedef struct AVD            AVD           ;
//typedef struct MarkerVolume   MarkerVolume  ;
//typedef struct PMatMono       PMatMono      ;
//typedef struct PMatBlock      PMatBlock     ;
//typedef struct MGLevel        MGLevel       ;
//typedef struct MG             MG            ;
//typedef struct WinStopCtx     WinStopCtx    ;




//---------------------------------------------------------------------------
#endif
