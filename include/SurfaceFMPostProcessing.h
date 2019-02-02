/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef SURFACEFMPOSTPROCESSING_H
#define SURFACEFMPOSTPROCESSING_H

// basic c++
#include <memory>

namespace sierra {
    namespace nalu {

        class PostProcessingData;
        class Realm;
        class SurfaceForceAndMomentAlgorithmDriver;

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

            std::unique_ptr<SurfaceForceAndMomentAlgorithmDriver> driverAlg_;
        };

    }  // nalu
}  // sierra


#endif /* SURFACEFMPOSTPROCESSING_H */
