import numpy as np
import time



from kimmdy.recipe import (
    Bind,
    Relax,
    Recipe,
    RecipeCollection,
)
from kimmdy.plugins import ReactionPlugin
from kimmdy.tasks import TaskFiles
from kimmdy.utils import (
    morse_transition_rate,
    get_atomnrs_from_plumedid,
    get_atominfo_from_atomnrs,
    get_bondprm_from_atomtypes,
    get_edissoc_from_atomnames,
)
from kimmdy.parsing import (
    read_top,
    read_plumed,
    read_distances_dat,
    read_edissoc,
)


class SpecBond(ReactionPlugin):
    """Bond formation. Implementation for time-varying rates within a cut-off
    """

    def get_recipe_collection(self, files: TaskFiles):
        logger = files.logger
        logger.debug("Getting recipe for reaction: specbond")

        # Initialization of filepaths
        files.input["itp"] = self.config.itp
        files.input["ebind"] = self.config.ebind
        frequency_factor = self.config.arrhenius_equation.frequency_factor
        dcutoff = self.config.dcutoff
        temperature = self.config.arrhenius_equation.temperature

        # Initialization of objects from files
        distances = read_distances_dat(files.input["plumed_out"])
        plumed = read_plumed(files.input["plumed"])
        top = self.runmng.top
        ffbonded = read_top(files.input["itp"])
        ebind = read_edissoc(files.input["ebind"]) #association energy (it uses the edissoc parsing functions )

        
        
        # handle already averaged distances
        t0, t1 = distances["time"][0], distances["time"][-1]
        
        if len(distances["time"]) == 1:
            logger.debug(
                "Plumed output contains single line, assuming averaged distance."
            )
            t0 = 0.0
            if t1 < 0.0001:
                logger.warning(
                    "First and last time in plumed output = 0.0"
                    "Check distance input.\n"
                    "If input is averaged already, set time to end time in ps."
                )
        if abs(t1 - t0) < 0.0001:
            logger.warning(
                "First and last time in plumed output equal!" "Check distance input."
            )

        recipes = []
        for plumedid, dists in distances.items():
            
            if plumedid == "time":
                continue
                
            #exclude pair if it is not within cutoff  ( dist[-1] = last time recorded by plumed )
            if float(dists[-1]) > dcutoff :
                continue

            
            atomnrs = get_atomnrs_from_plumedid(plumedid, plumed)
            
            #exclude bonds between adjacent beads    
            if  (int(atomnrs[0])-int(atomnrs[1]))**2 == 1 :
                continue


            
            #number of bonds
            Nbound=len( top.atoms[atomnrs[0]].bound_to_nrs )
            
            #decide if it will bind based on already formed bonds
            
            if (Nbound < 3):

                if(Nbound==1):  # first or last particle
                    attempt_bind=True  

                else:  #Nbonds = 2
                                     
                     idi=int ( atomnrs[0] )
                     idj1=int ( top.atoms[atomnrs[0]].bound_to_nrs[0] )
                     idj2=int (top.atoms[atomnrs[0]].bound_to_nrs[1] )
                     
                     
                     if ( (idj1 - idi)**2 == 1 ) and (  (idj2 - idi)**2 == 1 ) :  
                         attempt_bind=True
                         
                     else:
                         attempt_bind=False
                        
            else:  #Nbonds>=3
                attempt_bind=False
            
            
            

            # continue only if attempt bind = true 
            if attempt_bind==False:
                continue

                         
                           
            # get from plumedid to b0 and kb of the bond via atomtypes
            atomtypes, atomnames = get_atominfo_from_atomnrs(atomnrs, top)
            b0, kb = get_bondprm_from_atomtypes(atomtypes, ffbonded)
            residue = top.atoms[atomnrs[0]].residue
            Ebind = get_edissoc_from_atomnames(atomnames, ebind , residue)  #read association energy (Ebind for clarity)
                

            #Accept binding only with Boltzmann-weigthed probability
            R = 8.31446261815324e-3  # [kJ K-1 mol-1]
            Pbind=np.exp(-Ebind/(R * temperature ))
            k=frequency_factor*Pbind


                
            print( "Possible binding event: between atoms %s and %s with binding probability %f\n" %( atomnrs[0] , atomnrs[1], Pbind) )
                
                
            recipes.append(
                Recipe(
                    recipe_steps=[
                        Bind(atom_id_1=atomnrs[0], atom_id_2=atomnrs[1]),
                            
                    ],
                    rates=[(k)],
                    timespans=[(t0, t1)],
                )
            )
        
        return RecipeCollection(recipes)
