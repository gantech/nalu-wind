/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef SURFACEFMPOSTPROCESSING_H
#define SURFACEFMPOSTPROCESSING_H

#include <FieldTypeDef.h>

#include <stk_mesh/base/Part.hpp>

// basic c++
#include <string>
#include <memory>

namespace sierra {
    namespace nalu {

        class PostProcessingData;
        class Realm;

        struct SurfaceFMData {
            stk::mesh::PartVector partVector_;
            std::string outputFileName_;
            std::array<double,3> centroidCoords_;
            int frequency_;
            bool wallFunction;
        };
            
/** Post-processing to compute force and moment on various surfaces

 *  This class implements computing the force and moment on various surfaces

 *  Currently supported:
 *    - Compute pressure (pressureForce) and viscous force (tau_wall)
 *    - Use of wall function to compute viscous force in under-resolved 
 *      turbulent flow simulations
 *    - Computing yPlus based on distance of first node from the wall and 
 total tangential stress
*/
        
        class SurfaceFMPostProcessing
        {
        public:
            SurfaceFMPostProcessing(Realm&);

            ~SurfaceFMPostProcessing() = default;

            void register_surface_pp(const PostProcessingData&);

            void execute();

        private:
            Realm& realm_;

            std::vector<SurfaceFMData> surfaceFMData_;

            VectorFieldType *coordinates_;
            ScalarFieldType *pressure_;
            VectorFieldType *pressureForce_;
            ScalarFieldType *density_;
            ScalarFieldType *viscosity_;
            GenericFieldType *dudx_;
            VectorFieldType *tauWall_;
            ScalarFieldType *yplus_;
            GenericFieldType *exposedAreaVec_;
            ScalarFieldType *assembledArea_;
            ScalarFieldType *assembledAreaWF_;

            void zero_fields();

            void parallel_assemble_fields();
            
            void parallel_assemble_area();

            void cross_product(double *, double *, double *);

            void calc_surface_force(SurfaceFMData &);

            void create_file(std::string fileName);
        };

    }  // nalu
}  // sierra


#endif /* SURFACEFMPOSTPROCESSING_H */
