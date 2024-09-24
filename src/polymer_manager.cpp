
#include "polymer_manager.h"

const biosoup::NucleicAcid &compress_homopolymers(const biosoup::NucleicAcid &read) {
  std::string sequence = read.InflateData();
  std::string new_sequence = "";
  std::size_t seq_ptr = 0;
  seq_ptr += 1;
  std::size_t diamer_ptr = 0;
  diamer_ptr += 3;
  uint32_t repeat_count = 1;
  uint32_t direpeat_count = 1;
  std::vector<std::pair<std::size_t, std::size_t>> homopolymer_positions;
  std::vector<std::pair<std::size_t, std::size_t>> diamer_positions;
  std::size_t current_index = 0;
  std::string tmp = "";

  for(diamer_ptr; diamer_ptr<sequence.size();){
    tmp = sequence[diamer_ptr];
    if(sequence[diamer_ptr] == sequence[diamer_ptr-2] && sequence[diamer_ptr-1] == sequence[diamer_ptr-3] && !(sequence[diamer_ptr] == sequence[diamer_ptr-1])){
      repeat_count++;
      diamer_ptr+=2;
    }else{
      if(repeat_count > 1){
        diamer_positions.emplace_back(diamer_ptr-(2*repeat_count) - 1, diamer_ptr-2);
      }
      repeat_count = 1;
      diamer_ptr+=1;
    }
  }

  if(repeat_count > 1){
    diamer_positions.emplace_back(diamer_ptr-(2*repeat_count) - 1, diamer_ptr-2);
  }

  std::vector<std::pair<std::size_t, std::size_t>> new_diamer_pos;
  for(auto &positions : diamer_positions){
    new_sequence.append(sequence.substr(current_index, positions.first - current_index + 1));
    new_diamer_pos.emplace_back(new_sequence.size() - 1, new_sequence.size());
    current_index = positions.second;
  }
  if(current_index < sequence.size()){
    new_sequence.append(sequence.substr(current_index));
  }

  current_index = 0;
  sequence = "";
  std::size_t diamer_index = 0;
  repeat_count = 1;
  std::string tmp_2 = "";
  for(seq_ptr; seq_ptr < new_sequence.size();){
    tmp_2 = new_sequence[seq_ptr];
    if(new_sequence[seq_ptr] == new_sequence[seq_ptr-1] && (seq_ptr-1 != new_diamer_pos[diamer_index].second) && seq_ptr != new_diamer_pos[diamer_index].first){
      repeat_count++;
      seq_ptr++;
    }else if (seq_ptr == new_diamer_pos[diamer_index].first){
      if(diamer_index < new_diamer_pos.size()-1) diamer_index++;
      if(repeat_count > 1) homopolymer_positions.emplace_back(seq_ptr-repeat_count, seq_ptr-1);
      repeat_count = 1;
      seq_ptr += 3;

    }else if(seq_ptr == new_diamer_pos[diamer_index].second){
      if(diamer_index < new_diamer_pos.size()-1) diamer_index++;
      seq_ptr+=2;
      repeat_count = 1;
    }else {
      if(repeat_count > 1){
        homopolymer_positions.emplace_back(seq_ptr-repeat_count, seq_ptr-1);
      }
      repeat_count = 1;
      seq_ptr++;
    }
  }

  if(repeat_count > 1){
    homopolymer_positions.emplace_back(seq_ptr-repeat_count+1, seq_ptr-1);
  }

  std::vector<std::pair<std::size_t, std::size_t>> new_homopolymer_pos;
  std::vector<std::size_t> diamer_pos_final;
  current_index = 0;
  for(auto &positions : homopolymer_positions){
    sequence.append(new_sequence.substr(current_index, positions.first - current_index + 1));
    new_homopolymer_pos.emplace_back(sequence.size() - 1, sequence.size());
    current_index = positions.second + 1;
  }
  if(current_index < new_sequence.size()){
    sequence.append(new_sequence.substr(current_index));
  }

  return biosoup::NucleicAcid(read.name, sequence);
}