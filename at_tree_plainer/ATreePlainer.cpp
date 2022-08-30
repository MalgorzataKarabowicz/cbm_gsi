#include "ATreePlainer.hpp"

#include "AnalysisTree/TaskManager.hpp"

void ATreePlainer::Init()
{
  auto* man = AnalysisTree::TaskManager::GetInstance();
  auto* chain = man->GetChain();
  
  candidates_ = ANALYSISTREE_UTILS_GET<AnalysisTree::Particles*>(chain->GetPointerToBranch("Candidates"));
  simulated_= ANALYSISTREE_UTILS_GET<AnalysisTree::Particles*>(chain->GetPointerToBranch("Simulated"));
  can2sim_match_ = chain->GetMatchPointers().find(config_->GetMatchName("Candidates", "Simulated"))->second;
  
  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();
  
  AnalysisTree::BranchConfig out_particles("Complex", AnalysisTree::DetType::kParticle);
  out_particles.AddField<float>("chi2_geo");
  out_particles.AddField<float>("chi2_prim_first"); 
  out_particles.AddField<float>("chi2_prim_second");
  out_particles.AddField<float>("chi2_topo");       
  out_particles.AddField<float>("cosine_first");    
  out_particles.AddField<float>("cosine_second");   
  out_particles.AddField<float>("cosine_topo");     
  out_particles.AddField<float>("distance");        
  out_particles.AddField<float>("l");
  out_particles.AddField<float>("l_over_dl");
  out_particles.AddField<int>  ("generation");
  out_particles.AddField<int>  ("geant_process_id");
  out_particles.AddField<float>("lambda_chi2_geo");
  out_particles.AddField<float>("lambda_chi2_prim_first"); 
  out_particles.AddField<float>("lambda_chi2_prim_second");
  out_particles.AddField<float>("lambda_chi2_topo");       
  out_particles.AddField<float>("lambda_cosine_first");    
  out_particles.AddField<float>("lambda_cosine_second");   
  out_particles.AddField<float>("lambda_cosine_topo");     
  out_particles.AddField<float>("lambda_distance");        
  out_particles.AddField<float>("lambda_l");
  out_particles.AddField<float>("lambda_l_over_dl");
  out_particles.AddField<int>  ("lambda_generation");
  out_particles.AddField<int>  ("lambda_geant_process_id");
  out_particles.AddField<float>("lambda_mass_before_constraint");
  out_particles.AddField<float>("lambda_px");
  out_particles.AddField<float>("lambda_py");
  out_particles.AddField<float>("lambda_pz");
  
  man->AddBranch("Complex", plain_branch_, out_particles);
  
  if(cuts_)
    cuts_ -> Init(*out_config);
  
  InitIndices();
}

void ATreePlainer::Exec()
{
  plain_branch_ -> ClearChannels();
  
  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();
  
  for(auto& input_particle : *candidates_)
  {
    if(cuts_)
      if(!cuts_->Apply(input_particle))
        continue;
            
    auto& output_particle = plain_branch_->AddChannel(out_config->GetBranchConfig(plain_branch_->GetId()));
    output_particle.SetMomentum(input_particle.GetPx(), input_particle.GetPy(), input_particle.GetPz());
    output_particle.SetMass(input_particle.GetMass());
    output_particle.SetPid(input_particle.GetPid());
    output_particle.SetField(input_particle.GetField<float>(chi2_geo_field_id_r_        ), chi2_geo_field_id_w1_        );
    output_particle.SetField(input_particle.GetField<float>(chi2_prim_first_field_id_r_ ), chi2_prim_first_field_id_w1_ );
    output_particle.SetField(input_particle.GetField<float>(chi2_prim_second_field_id_r_), chi2_prim_second_field_id_w1_);
    output_particle.SetField(input_particle.GetField<float>(chi2_topo_field_id_r_       ), chi2_topo_field_id_w1_       );
    output_particle.SetField(input_particle.GetField<float>(cosine_first_field_id_r_    ), cosine_first_field_id_w1_    );
    output_particle.SetField(input_particle.GetField<float>(cosine_second_field_id_r_   ), cosine_second_field_id_w1_   );
    output_particle.SetField(input_particle.GetField<float>(cosine_topo_field_id_r_     ), cosine_topo_field_id_w1_     );
    output_particle.SetField(input_particle.GetField<float>(distance_field_id_r_        ), distance_field_id_w1_        );
    output_particle.SetField(input_particle.GetField<float>(l_field_id_r_               ), l_field_id_w1_               );
    output_particle.SetField(input_particle.GetField<float>(l_over_dl_field_id_r_       ), l_over_dl_field_id_w1_       );
    output_particle.SetField(input_particle.GetField<int  >(generation_field_id_r_      ), generation_field_id_w1_      );
    
    if(input_particle.GetField<int  >(generation_field_id_r_)>0)
    {
      auto& matched_particle = simulated_->GetChannel(can2sim_match_->GetMatch(input_particle.GetId()));
      output_particle.SetField(matched_particle.GetField<int>(g4_process_id_field_id_r_), g4_process_id_field_id_w1_);
    }
    else
      output_particle.SetField(-999, g4_process_id_field_id_w1_);
        
    //-------------------------------------------------------------------------------------------------------------------------------
    const int lambda_id = input_particle.GetField<int>(daughter_id2_field_id_r_);
    auto& additional_particle = candidates_->GetChannel(lambda_id);
    
    output_particle.SetField(additional_particle.GetField<float>(chi2_geo_field_id_r_        ), chi2_geo_field_id_w2_        );
    output_particle.SetField(additional_particle.GetField<float>(chi2_prim_first_field_id_r_ ), chi2_prim_first_field_id_w2_ );
    output_particle.SetField(additional_particle.GetField<float>(chi2_prim_second_field_id_r_), chi2_prim_second_field_id_w2_);
    output_particle.SetField(additional_particle.GetField<float>(chi2_topo_field_id_r_       ), chi2_topo_field_id_w2_       );
    output_particle.SetField(additional_particle.GetField<float>(cosine_first_field_id_r_    ), cosine_first_field_id_w2_    );
    output_particle.SetField(additional_particle.GetField<float>(cosine_second_field_id_r_   ), cosine_second_field_id_w2_   );
    output_particle.SetField(additional_particle.GetField<float>(cosine_topo_field_id_r_     ), cosine_topo_field_id_w2_     );
    output_particle.SetField(additional_particle.GetField<float>(distance_field_id_r_        ), distance_field_id_w2_        );
    output_particle.SetField(additional_particle.GetField<float>(l_field_id_r_               ), l_field_id_w2_               );
    output_particle.SetField(additional_particle.GetField<float>(l_over_dl_field_id_r_       ), l_over_dl_field_id_w2_       );
    output_particle.SetField(additional_particle.GetField<int  >(generation_field_id_r_      ), generation_field_id_w2_      );
    output_particle.SetField(additional_particle.GetPx(), px_field_id_w2_);
    output_particle.SetField(additional_particle.GetPy(), py_field_id_w2_);
    output_particle.SetField(additional_particle.GetPz(), pz_field_id_w2_);
    
    const float mass_before_constraint = lambda_mass + additional_particle.GetField<float>(invmass_discr_field_id_r_)*lambda_mass_sigma;
    output_particle.SetField(mass_before_constraint, mass_before_constraint_field_id_w2_);
        
    if(additional_particle.GetField<int  >(generation_field_id_r_)>0)
    {
      auto& matched_particle = simulated_->GetChannel(can2sim_match_->GetMatch(additional_particle.GetId()));
      output_particle.SetField(matched_particle.GetField<int>(g4_process_id_field_id_r_), g4_process_id_field_id_w2_);
    }
    else
      output_particle.SetField(-999, g4_process_id_field_id_w2_);    
  }
  
}

void ATreePlainer::InitIndices()
{
  auto in_branch_cand = config_->GetBranchConfig("Candidates");
  auto in_branch_sim  = config_->GetBranchConfig("Simulated");
  
  chi2_geo_field_id_r_         = in_branch_cand.GetFieldId("chi2_geo");
  chi2_prim_first_field_id_r_  = in_branch_cand.GetFieldId("chi2_prim_first"); 
  chi2_prim_second_field_id_r_ = in_branch_cand.GetFieldId("chi2_prim_second");
  chi2_topo_field_id_r_        = in_branch_cand.GetFieldId("chi2_topo");       
  cosine_first_field_id_r_     = in_branch_cand.GetFieldId("cosine_first");    
  cosine_second_field_id_r_    = in_branch_cand.GetFieldId("cosine_second");   
  cosine_topo_field_id_r_      = in_branch_cand.GetFieldId("cosine_topo");     
  distance_field_id_r_         = in_branch_cand.GetFieldId("distance");        
  l_field_id_r_                = in_branch_cand.GetFieldId("l");
  l_over_dl_field_id_r_        = in_branch_cand.GetFieldId("l_over_dl");
  daughter_id1_field_id_r_     = in_branch_cand.GetFieldId("daughter1_id");
  daughter_id2_field_id_r_     = in_branch_cand.GetFieldId("daughter2_id");
  generation_field_id_r_       = in_branch_cand.GetFieldId("generation");
  invmass_discr_field_id_r_    = in_branch_cand.GetFieldId("invmass_discr");
  
  g4_process_id_field_id_r_    = in_branch_sim.GetFieldId("geant_process_id");
  
  
  auto out_config = AnalysisTree::TaskManager::GetInstance()->GetConfig();
  const auto& out_branch = out_config->GetBranchConfig(plain_branch_->GetId());
  
  chi2_geo_field_id_w1_               = out_branch.GetFieldId("chi2_geo");
  chi2_prim_first_field_id_w1_        = out_branch.GetFieldId("chi2_prim_first"); 
  chi2_prim_second_field_id_w1_       = out_branch.GetFieldId("chi2_prim_second");
  chi2_topo_field_id_w1_              = out_branch.GetFieldId("chi2_topo");       
  cosine_first_field_id_w1_           = out_branch.GetFieldId("cosine_first");    
  cosine_second_field_id_w1_          = out_branch.GetFieldId("cosine_second");   
  cosine_topo_field_id_w1_            = out_branch.GetFieldId("cosine_topo");     
  distance_field_id_w1_               = out_branch.GetFieldId("distance");        
  l_field_id_w1_                      = out_branch.GetFieldId("l");
  l_over_dl_field_id_w1_              = out_branch.GetFieldId("l_over_dl");
  generation_field_id_w1_             = out_branch.GetFieldId("generation");
  g4_process_id_field_id_w1_          = out_branch.GetFieldId("geant_process_id");
  chi2_geo_field_id_w2_               = out_branch.GetFieldId("lambda_chi2_geo");
  chi2_prim_first_field_id_w2_        = out_branch.GetFieldId("lambda_chi2_prim_first"); 
  chi2_prim_second_field_id_w2_       = out_branch.GetFieldId("lambda_chi2_prim_second");
  chi2_topo_field_id_w2_              = out_branch.GetFieldId("lambda_chi2_topo");       
  cosine_first_field_id_w2_           = out_branch.GetFieldId("lambda_cosine_first");    
  cosine_second_field_id_w2_          = out_branch.GetFieldId("lambda_cosine_second");   
  cosine_topo_field_id_w2_            = out_branch.GetFieldId("lambda_cosine_topo");     
  distance_field_id_w2_               = out_branch.GetFieldId("lambda_distance");        
  l_field_id_w2_                      = out_branch.GetFieldId("lambda_l");
  l_over_dl_field_id_w2_              = out_branch.GetFieldId("lambda_l_over_dl");
  generation_field_id_w2_             = out_branch.GetFieldId("lambda_generation");
  g4_process_id_field_id_w2_          = out_branch.GetFieldId("lambda_geant_process_id");
  mass_before_constraint_field_id_w2_ = out_branch.GetFieldId("lambda_mass_before_constraint");
  px_field_id_w2_                     = out_branch.GetFieldId("lambda_px");
  py_field_id_w2_                     = out_branch.GetFieldId("lambda_py");
  pz_field_id_w2_                     = out_branch.GetFieldId("lambda_pz");
}