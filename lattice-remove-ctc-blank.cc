// MIT License
//
// Copyright (c) 2016 Joan Puigcerver <joapuipe@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "base/kaldi-common.h"
#include "util/common-utils.h"
#include "fstext/kaldi-fst-io.h"
#include "fstext/fstext-utils.h"
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"

namespace kaldi {

void RemoveCTCBlankFromLattice(
    const Lattice& inp, const LatticeArc::Label blank, Lattice* out) {
  typedef LatticeArc::StateId StateId;
  typedef LatticeArc::Label Label;
  // Get table mapping from output symbols from the output labels in the input
  // lattice to states
  std::unordered_map<Label, StateId> symbol2state;
  for (fst::StateIterator<Lattice> siter(inp); !siter.Done(); siter.Next()) {
    const StateId s = siter.Value();
    for (fst::ArcIterator<Lattice> aiter(inp, s); !aiter.Done();
         aiter.Next()) {
      const LatticeArc& arc = aiter.Value();
      const Label o = arc.olabel;
      if (o != blank && o != 0 && symbol2state.count(o) == 0) {
        symbol2state.insert(std::make_pair(o, symbol2state.size() + 1));
      }
    }
  }
  // Create composition lattice, such that output = compose(input, C)
  Lattice C;
  for (size_t s = 0; s < symbol2state.size() + 1; ++s) {
    C.SetFinal(C.AddState(), LatticeWeight::One());
  }
  C.SetStart(0);
  // Self-loop in the blank state
  C.AddArc(0, LatticeArc(blank, 0, LatticeWeight::One(), 0));
  for (const std::pair<Label, StateId>& p : symbol2state) {
    // Arc from the initial state to the symbol's state
    C.AddArc(0, LatticeArc(p.first, p.first, LatticeWeight::One(), p.second));
    // Self-loop in the symbol's state
    C.AddArc(p.second, LatticeArc(p.first, 0, LatticeWeight::One(), p.second));
    // Arc back to the blank state
    C.AddArc(p.second, LatticeArc(blank, 0, LatticeWeight::One(), 0));
    // Arc to all other symbol states
    for (const std::pair<Label, StateId>& p2 : symbol2state) {
      if (p.first != p2.first) {
        C.AddArc(p.second, LatticeArc(p2.first, p2.first,
                                      LatticeWeight::One(), p2.second));
      }
    }
  }
  // Compute output lattice
  fst::Compose(inp, C, out);
}

}  // namespace fst

int main(int argc, char** argv) {
  try {
    using namespace kaldi;

    const char* usage =
        "Remove CTC blank symbols from the output labels of Kaldi lattices.\n"
        "\n"
        "Usage: lattice-remove-ctc-blank blank-symbol lat-rspecifier lat-wspecifier\n"
        " e.g.: lattice-remove-ctc-blank 32 ark:input.ark ark:output.ark\n";

    ParseOptions po(usage);
    po.Read(argc, argv);

    if (po.NumArgs() != 3) {
      po.PrintUsage();
      exit(1);
    }

    const std::string blank_symbol_str = po.GetArg(1);
    const std::string lattice_in_str = po.GetArg(2);
    const std::string lattice_out_str = po.GetArg(3);
    const bool lattice_in_is_table =
        (ClassifyRspecifier(lattice_in_str, NULL, NULL) != kNoRspecifier);
    const bool lattice_out_is_table =
        (ClassifyRspecifier(lattice_out_str, NULL, NULL) != kNoRspecifier);

    LatticeArc::Label blank_symbol = 0;
    if (!ConvertStringToInteger(blank_symbol_str, &blank_symbol)) {
      KALDI_ERR << "String \"" << blank_symbol_str
                << "\" cannot be converted to an integer";
    }
    if (blank_symbol == 0) {
      KALDI_ERR << "Symbol 0 is reserved for epsilon!";
    }


    if (lattice_in_is_table && lattice_out_is_table) {
      SequentialLatticeReader lattice_reader(lattice_in_str);
      LatticeWriter lattice_writer(lattice_out_str);
      for (; !lattice_reader.Done(); lattice_reader.Next()) {
        // Read input lattice
        const std::string lattice_key = lattice_reader.Key();
        Lattice lat = lattice_reader.Value();
        lattice_reader.FreeCurrent();
        // Make sure that lattice complies with all asumptions
        const uint64_t properties =
            lat.Properties(fst::kAcceptor | fst::kAcyclic, true);
        if ((properties & fst::kAcceptor) != fst::kAcceptor) {
          KALDI_ERR << "Lattice " << lattice_key << " is not an acceptor";
        }
        if ((properties & fst::kAcyclic) != fst::kAcyclic) {
          KALDI_ERR << "Lattice " << lattice_key << " is not acyclic";
        }
        Lattice out;
        RemoveCTCBlankFromLattice(lat, blank_symbol, &out);
        lattice_writer.Write(lattice_key, out);
      }
    } else {
      KALDI_ERR << "Not implemented! Both input and output lattices must be "
                << "Kaldi tables.";
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what();
    return 1;
  }
}
