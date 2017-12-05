#include <iostream>
#include <map>
#include <sbml/SBMLTypes.h>

using namespace std;
LIBSBML_CPP_NAMESPACE_USE

ASTNode* addASTasReactant(ASTNode* ast, Reaction* r) {
	ASTNode *root, *l;
	if (ast == nullptr) {  // If there is no parent, return -1 * v1.
		root = new ASTNode(AST_TIMES);
		l = new ASTNode();
		l->setValue(-1.0);
		root->addChild(l);
		root->addChild(ASTNode_deepCopy(r->getKineticLaw()->getMath()));
	} else {
		root = new ASTNode(AST_MINUS);
		root->addChild(ast);
		root->addChild(ASTNode_deepCopy(r->getKineticLaw()->getMath()));
	}
	return root;
}

ASTNode* addASTasProduct(ASTNode* ast, Reaction* r) {
	ASTNode *root;
	if (ast == nullptr) {  // If there is no parent, just return v1.
		root = ASTNode_deepCopy(r->getKineticLaw()->getMath());
	} else {
		root = new ASTNode(AST_PLUS);
		root->addChild(ast);
		root->addChild(ASTNode_deepCopy(r->getKineticLaw()->getMath()));
	}
	return root;
}

bool isSpeciesReactantOf(Species* s, Reaction* r) {
	for (auto i = 0; i < r->getNumReactants(); i++) {
		auto reactant = r->getReactant(i);
		if (reactant->getSpecies() == s->getId()) {
			return true;
		}
	}
	return false;
}

bool isSpeciesProductOf(Species* s, Reaction* r) {
	for (auto i = 0; i < r->getNumProducts(); i++) {
		auto product = r->getProduct(i);
		if (product->getSpecies() == s->getId()) {
			return true;
		}
	}
	return false;
}

map<string, ASTNode*> loadModel(const char* filename) {
	SBMLReader reader;

	auto document = reader.readSBML(filename);
	auto model = document->getModel();

	map<string, ASTNode*> hashRate;

	// Generate Rate equation for all variable Species (ex. dx/dt = v1 - v2 + v3).
	for (auto i = 0; i < model->getNumSpecies(); i++) {
		auto species = model->getSpecies(i);
		if (species->getBoundaryCondition() || species->getConstant()) {
			continue;
		}
		ASTNode* root = nullptr;
		for (auto j = 0; j < model->getNumReactions(); j++) {
			auto reaction = model->getReaction(j);
			if (isSpeciesReactantOf(species, reaction)) {
				root = addASTasReactant(root, reaction);
			}
			if (isSpeciesProductOf(species, reaction)) {
				root = addASTasProduct(root, reaction);
			}
		}
		if (root != nullptr) {
			hashRate[species->getId()] = root;
		}
	}
  return hashRate;
}

void printRHS(map<string, ASTNode*> hashRate) {
	for(auto itr = hashRate.begin(); itr != hashRate.end(); ++itr) {
		cout << "d[" << itr->first << "]/dt = " << SBML_formulaToString(itr->second) << endl;
	}
}

int main(int argc, char const* argv[])
{
  if (argc != 2)
  {
    cout << endl << "Usage: " << argv[0] << " filename" << endl << endl;
    return 1;
  }
  const char* filename  = argv[1];
  SBMLReader reader;

  auto document = reader.readSBML(filename);
	auto hashRate = loadModel(filename);

  printRHS(hashRate);

	delete document;
	return 0;
}
