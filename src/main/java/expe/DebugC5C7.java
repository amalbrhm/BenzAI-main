package expe;// Dans votre classe de débogage, par exemple src/main/java/expe/DebugC5C7.java

import generator.GeneralModel;
import generator.ModelBuilder;
import generator.SolverResults;
import generator.properties.model.ModelPropertySet;
import generator.properties.model.expression.BinaryNumericalExpression;

public class DebugC5C7 {

    public static void main(String[] args) {
        System.out.println("--- Début Débogage CSP pour C5/C7 (Étape Initiale) ---");

        ModelPropertySet modelPropertySet = new ModelPropertySet();

        // Configuration minimale pour tester :
        // Forcer une petite grille, et demander 0 C5, 0 C7 pour commencer.
        // Cela devrait se comporter comme une génération normale d'hexagones C6.
        //modelPropertySet.getById("coronenoid").addExpression(new BinaryNumericalExpression("coronenoid", "=", 2));
        modelPropertySet.getById("hexagons").addExpression(new BinaryNumericalExpression("hexagons", "=", 7));
        modelPropertySet.getById("pentagons").addExpression(new BinaryNumericalExpression("pentagons", "=", 2));
        modelPropertySet.getById("heptagons").addExpression(new BinaryNumericalExpression("heptagons", "=", 2));

        // Optionnel : forcez un nombre de couronnes, par ex. 2 couronnes pour avoir quelques sites potentiels.
        // (2 couronnes = 7 hexagones de base au total dans la grille)
        // Si "hexagons" = 2 et transformations = 0, il faut une grille de base d'au moins 2 hexagones.
        // Si on veut des transformations, il faut une grille de base plus grande.
        // Ajustez ceci pour avoir un nombre raisonnable de `potentialTransformationSites`.
        // Par exemple, pour avoir des sites, il faut au moins 2 hexagones de base adjacents.
        // Si vous demandez peu d'hexagones finaux (ex: 2 C6) et 0 transformations,
        // une petite grille de base (ex: 2 couronnes) est bien.
        // Forcer explicitement le nombre de couronnes pour le test :


        GeneralModel model = ModelBuilder.buildModel(modelPropertySet);
        if (model == null) {
            System.err.println("Le modèle GeneralModel n'a pas pu être construit.");
            return;
        }
        model.setInTestMode(true); // TRES IMPORTANT POUR LE DEBOGAGE

        // Afficher les infos du modèle AVANT résolution
        System.out.println("\nConfiguration du GeneralModel AVANT résolution:");
        System.out.println("  - Nb Couronnes grille de base: " + model.getNbCrowns());
        System.out.println("  - Nb Hexagones grille de base (nbHexagonsCoronenoid): " + model.getNbHexagonsCoronenoid());
        if (model.getPotentialTransformationSites() != null) {
            System.out.println("  - Nb Sites de Transformation Potentiels: " + model.getPotentialTransformationSites().size());
        } else {
            System.out.println("  - Nb Sites de Transformation Potentiels: null (problème d'initialisation ?)");
        }
        if (model.getTransformationActiveVars() != null) {
            System.out.println("  - Nb Variables transformationActiveVars: " + model.getTransformationActiveVars().length);
        } else {
            System.out.println("  - transformationActiveVars est null (problème d'initialisation ?)");
        }
        if (model.getNbPentagonsVar() != null) { // Doit être initialisé par PentagonNumberConstraint.buildVariables()
            System.out.println("  - Domaine initial nbPentagonsVar: [" + model.getNbPentagonsVar().getLB() + "," + model.getNbPentagonsVar().getUB() + "]");
        } else {
            System.out.println("  - nbPentagonsVar est null dans GeneralModel (non défini par contrainte ?)");
        }
        if (model.getNbHeptagonsVar() != null) {
            System.out.println("  - Domaine initial nbHeptagonsVar: [" + model.getNbHeptagonsVar().getLB() + "," + model.getNbHeptagonsVar().getUB() + "]");
        } else {
            System.out.println("  - nbHeptagonsVar est null dans GeneralModel (non défini par contrainte ?)");
        }

        System.out.println("\nLancement de la résolution Choco...");
        SolverResults results = model.solve(); // La méthode solve() contient displaySolution()

        System.out.println("\n--- Fin de la Résolution Choco ---");
        System.out.println("Nombre total de solutions Choco (configurations) trouvées: " + results.getNbTotalSolutions());
        System.out.println("Temps de résolution: " + results.getTime() + "ms");

        // La méthode displaySolution() dans GeneralModel devrait déjà afficher les valeurs des variables pour chaque solution.
    }
}