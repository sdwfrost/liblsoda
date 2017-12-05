#include <iostream>
#include <map>
#include <sbml/SBMLTypes.h>

extern "C" {
#include "lsoda.h"
}

using namespace std;
LIBSBML_CPP_NAMESPACE_USE

typedef struct data {
  map<string, ASTNode *> mapRates;
  map<string, double> mapVariables;
  Model *model;
} data_t;

ASTNode *addASTasReactant(ASTNode *ast, ASTNode *kineticlaw) {
  ASTNode *root, *l;
  if (ast == nullptr) {  // If there is no parent, return -1 * v1.
    root = new ASTNode(AST_TIMES);
    l = new ASTNode();
    l->setValue(-1.0);
    root->addChild(l);
    root->addChild(kineticlaw);
  } else {
    root = new ASTNode(AST_MINUS);
    root->addChild(ast);
    root->addChild(kineticlaw);
  }
  return root;
}

ASTNode *addASTasProduct(ASTNode *ast, ASTNode *kineticlaw) {
  ASTNode *root;
  if (ast == nullptr) {  // If there is no parent, just return v1.
    root = kineticlaw;
  } else {
    root = new ASTNode(AST_PLUS);
    root->addChild(ast);
    root->addChild(kineticlaw);
  }
  return root;
}

bool isSpeciesReactantOf(Species *s, Reaction *r) {
  for (auto i = 0; i < r->getNumReactants(); i++) {
    auto reactant = r->getReactant(i);
    if (reactant->getSpecies() == s->getId()) {
      return true;
    }
  }
  return false;
}

bool isSpeciesProductOf(Species *s, Reaction *r) {
  for (auto i = 0; i < r->getNumProducts(); i++) {
    auto product = r->getProduct(i);
    if (product->getSpecies() == s->getId()) {
      return true;
    }
  }
  return false;
}

ASTNode *rewriteLocalParameters(const ASTNode *node, const ListOfParameters *localParameters) {
  ASTNode *ret;

  if (node->getType() == AST_NAME) {
    auto name = node->getName();
    for (auto i = 0; i < localParameters->size(); i++) {
      auto param = localParameters->get(i);
      if (name == param->getId()) {
        ret = new ASTNode(AST_REAL);
        ret->setValue(param->getValue());
        return ret;
      }
    }
  }

  // if the node doesn't represent local parameter

  ret = node->deepCopy();

  // replace children's local parameter node recursively
  for (auto i = 0; i < ret->getNumChildren(); i++) {
    auto newChild = rewriteLocalParameters(ret->getChild(i), localParameters);
    ret->replaceChild(i, newChild, true);
  }

  return ret;
}

data_t loadModel(const char *filename) {
  data_t mydata;
  SBMLReader reader;

  auto document = reader.readSBML(filename);
  mydata.model = document->getModel();

  // rewrite local parameters
  map<string, ASTNode *> mapKineticLaw;
  for (auto i = 0; i < mydata.model->getNumReactions(); i++) {
    auto reaction = mydata.model->getReaction(i);
    auto node = reaction->getKineticLaw()->getMath();
    mapKineticLaw[reaction->getId()] = rewriteLocalParameters(node, reaction->getKineticLaw()->getListOfParameters());
  }

  // Generate Rate equation and map for all variable Species (ex. dx/dt = v1 - v2 + v3).
  for (auto i = 0; i < mydata.model->getNumSpecies(); i++) {
    auto species = mydata.model->getSpecies(i);

    // if constant
    if (species->getBoundaryCondition() || species->getConstant()) {
      continue;
    }

    // InitialAmount or InitialConcentration
    if (species->isSetInitialAmount()) {
      mydata.mapVariables[species->getId()] = species->getInitialAmount();
    } else {
      mydata.mapVariables[species->getId()] = species->getInitialConcentration();
    }

    // generate Rate equation
    ASTNode *root = nullptr;
    for (auto j = 0; j < mydata.model->getNumReactions(); j++) {
      auto reaction = mydata.model->getReaction(j);
      if (isSpeciesReactantOf(species, reaction)) {
        root = addASTasReactant(root, mapKineticLaw[reaction->getId()]);
      }
      if (isSpeciesProductOf(species, reaction)) {
        root = addASTasProduct(root, mapKineticLaw[reaction->getId()]);
      }
    }
    if (root != nullptr) {
      mydata.mapRates[species->getId()] = root;
    }
  }
  SBMLDocument_free(document);
  return mydata;
}

void printRHS(map<string, ASTNode *> mapRate) {
  cout << "[Rate equations]" << endl;
  for (auto itr : mapRate) {
    cout << "d[" << itr.first << "]/dt = " << SBML_formulaToString(itr.second) << endl;
  }
  cout << endl;
}

void printVariables(map<string, double> mapVariables) {
  cout << "[Species]" << endl;
  for (auto itr : mapVariables) {
    cout << itr.first << " = " << itr.second << endl;
  }
  cout << endl;
}

int fex(double t, double *y, double *ydot, void *data) {
  data_t *data_p = (data_t *) data;
  int idx = 0;
  for (auto itr : data_p->mapRates) {
    string key = itr.first;
    ASTNode *ast = itr.second;
    ydot[idx] = SBMLTransforms::evaluateASTNode(ast, data_p->mapVariables, data_p->model);
    idx++;
  }
  return 0;
}

void printResult(double t, double *y, struct lsoda_context_t* ctx) {
  cout << t << "\t";
  for (auto i = 0; i < ctx->neq; i++) {
    cout << y[i] << "\t";
  }
  cout << endl;

}

void printResultHeader(double t, double *y, struct lsoda_context_t* ctx) {
  cout << "time" << "\t";
  for (auto itr : ((data_t*)ctx->data)->mapVariables) {
    cout << itr.first << "\t";
  }
  cout << endl;
  printResult(t, y, ctx);
}

int main(int argc, char const *argv[]) {
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " filename" << endl;
    return 1;
  }
  const char *filename = argv[1];
  auto mydata = loadModel(filename);
  printVariables(mydata.mapVariables);
  printRHS(mydata.mapRates);

  // LSODA
  int neq = mydata.mapRates.size();
  double *rtol = new double[neq];
  double *atol = new double[neq];
  double *y = new double[neq];

  for (auto i = 0; i < neq; i++) {
    rtol[i] = 1e-6;
    atol[i] = 1e-6;
  }
  int idx = 0;
  for (auto itr : mydata.mapVariables) {
    y[idx] = itr.second;
    idx++;
  }

  struct lsoda_opt_t opt = {0};
  opt.ixpr = 0;
  opt.rtol = rtol;
  opt.atol = atol;
  opt.itask = 1;

  struct lsoda_context_t ctx = {
      .function = fex,
      .neq = neq,
      .data = &mydata,
      .state = 1,
  };

  lsoda_prepare(&ctx, &opt);

  double t = 0.0;
  printResultHeader(t, y, &ctx);
  for (double tout = 0.1; tout <= 1.0; tout += 0.1) {
    lsoda(&ctx, y, &t, tout);
    printResult(t, y, &ctx);

    if (ctx.state <= 0) {
      std::cout << "[ERROR] ctx.stat <= 0" << std::endl;
    }
  }

  lsoda_free(&ctx);

  return 0;
}
