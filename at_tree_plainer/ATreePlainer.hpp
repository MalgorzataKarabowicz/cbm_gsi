#ifndef ATREEPLAINER_HPP
#define ATREEPLAINER_HPP

#include "AnalysisTree/Cuts.hpp"
#include "AnalysisTree/Detector.hpp"
#include "AnalysisTree/Task.hpp"

class ATreePlainer : public AnalysisTree::Task {
 public:
  explicit ATreePlainer() = default;
  ~ATreePlainer() override = default;

  void Init() override;
  void Exec() override;
  void Finish() override {};
  
  void InitIndices();
  
  void SetCuts(AnalysisTree::Cuts* cuts) { cuts_ = cuts; }
  
  
protected:
  
  // input branches
  AnalysisTree::Particles* candidates_{nullptr};
  AnalysisTree::Particles* simulated_{nullptr};
  AnalysisTree::Matching* can2sim_match_{nullptr};
  
  // output branch
  AnalysisTree::Particles* plain_branch_{nullptr};
  
  AnalysisTree::Cuts* cuts_{nullptr};
  
  static constexpr float lambda_mass = 1.115683;
  static constexpr float lambda_mass_sigma = 1.5e-3;
  
  //**** input fields ***********
  int chi2_geo_field_id_r_{AnalysisTree::UndefValueInt};
  int chi2_prim_first_field_id_r_{AnalysisTree::UndefValueInt};  
  int chi2_prim_second_field_id_r_{AnalysisTree::UndefValueInt}; 
  int chi2_topo_field_id_r_{AnalysisTree::UndefValueInt};        
  int cosine_first_field_id_r_{AnalysisTree::UndefValueInt};     
  int cosine_second_field_id_r_{AnalysisTree::UndefValueInt};    
  int cosine_topo_field_id_r_{AnalysisTree::UndefValueInt};      
  int distance_field_id_r_{AnalysisTree::UndefValueInt};        
  int l_field_id_r_{AnalysisTree::UndefValueInt};
  int l_over_dl_field_id_r_{AnalysisTree::UndefValueInt};
  int daughter_id1_field_id_r_{AnalysisTree::UndefValueInt};
  int daughter_id2_field_id_r_{AnalysisTree::UndefValueInt};
  int generation_field_id_r_{AnalysisTree::UndefValueInt};
  int invmass_discr_field_id_r_{AnalysisTree::UndefValueInt};  
  
  int g4_process_id_field_id_r_{AnalysisTree::UndefValueInt};
  //*****************************
  
  //***** output fields *********
  int chi2_geo_field_id_w1_{AnalysisTree::UndefValueInt};
  int chi2_prim_first_field_id_w1_{AnalysisTree::UndefValueInt};  
  int chi2_prim_second_field_id_w1_{AnalysisTree::UndefValueInt}; 
  int chi2_topo_field_id_w1_{AnalysisTree::UndefValueInt};        
  int cosine_first_field_id_w1_{AnalysisTree::UndefValueInt};     
  int cosine_second_field_id_w1_{AnalysisTree::UndefValueInt};    
  int cosine_topo_field_id_w1_{AnalysisTree::UndefValueInt};      
  int distance_field_id_w1_{AnalysisTree::UndefValueInt};        
  int l_field_id_w1_{AnalysisTree::UndefValueInt};
  int l_over_dl_field_id_w1_{AnalysisTree::UndefValueInt};
  int generation_field_id_w1_{AnalysisTree::UndefValueInt};
  int g4_process_id_field_id_w1_{AnalysisTree::UndefValueInt};
  //-----------------------------
  int chi2_geo_field_id_w2_{AnalysisTree::UndefValueInt};
  int chi2_prim_first_field_id_w2_{AnalysisTree::UndefValueInt};  
  int chi2_prim_second_field_id_w2_{AnalysisTree::UndefValueInt}; 
  int chi2_topo_field_id_w2_{AnalysisTree::UndefValueInt};        
  int cosine_first_field_id_w2_{AnalysisTree::UndefValueInt};     
  int cosine_second_field_id_w2_{AnalysisTree::UndefValueInt};    
  int cosine_topo_field_id_w2_{AnalysisTree::UndefValueInt};      
  int distance_field_id_w2_{AnalysisTree::UndefValueInt};        
  int l_field_id_w2_{AnalysisTree::UndefValueInt};
  int l_over_dl_field_id_w2_{AnalysisTree::UndefValueInt};
  int generation_field_id_w2_{AnalysisTree::UndefValueInt};
  int g4_process_id_field_id_w2_{AnalysisTree::UndefValueInt};
  int mass_before_constraint_field_id_w2_{AnalysisTree::UndefValueInt};
  int px_field_id_w2_{AnalysisTree::UndefValueInt};
  int py_field_id_w2_{AnalysisTree::UndefValueInt};
  int pz_field_id_w2_{AnalysisTree::UndefValueInt};
  //******************************
  
};
#endif // ATREEPLAINER_HPP