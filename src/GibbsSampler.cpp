//
// Created by pietr on 12/12/2021.
//

#include "GibbsSampler.h"


void GibbsSampler::GS_Step() {
  for(FullConditional* full_cond: this->FullConditionals){ // mettere prima update della partition (da aggiungere anche prima)
    std::cout<<1<<std::endl;
    //full_cond->update(gs_data, random_engine);
    std::cout<<full_cond->name<<std::endl;
  }
  //std::cout << gs_data.M << "\n";
  //std::cout << gs_data.K << "\n";
}

void GibbsSampler::store_params_values() {
    for(std::map<string,std::vector<double>>::iterator iter = output_data.begin(); iter != output_data.end(); ++iter)
        output_data[iter->first].push_back(parameters[iter->first]);
}


std::map<string, std::vector<double>> GibbsSampler::sample() {
    for(unsigned int it; it<burn_in + n_iter * thin; it++){
        this->GS_Step();
        if(it>burn_in && it%thin == 0){
            parameters["M"]=gs_data.M;
            parameters["M*"]=gs_data.Mstar;
            parameters["K"]=gs_data.K;
            this->store_params_values();
        }
    }
    return output_data;
}
