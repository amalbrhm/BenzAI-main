package generator;

import generator.patterns.Pattern;
import generator.patterns.PatternLabel;
import generator.patterns.PatternOccurences;
import generator.properties.Property;
import generator.properties.model.ModelProperty;
import generator.properties.model.ModelPropertySet;
import generator.properties.model.expression.ParameterizedExpression;
import generator.properties.solver.SolverProperty;
import generator.properties.solver.SolverPropertySet;
import javafx.application.Platform;
import javafx.beans.property.SimpleIntegerProperty;
import benzenoid.Benzenoid;
import benzenoid.Node;
import nogood.*;
import org.chocosolver.solver.Model;
import org.chocosolver.solver.Solver;
import org.chocosolver.solver.search.strategy.selectors.values.IntDomainMax;
import org.chocosolver.solver.search.strategy.selectors.values.IntDomainMin;
import org.chocosolver.solver.search.strategy.selectors.variables.FirstFail;
import org.chocosolver.solver.search.strategy.strategy.IntStrategy;
import org.chocosolver.solver.variables.BoolVar;
import org.chocosolver.solver.variables.IntVar;
import org.chocosolver.solver.variables.UndirectedGraphVar;
import org.chocosolver.solver.variables.Variable;
import org.chocosolver.util.objects.graphs.UndirectedGraph;
import solution.BenzenoidSolution;
import utils.Couple;
import utils.HexNeighborhood;
import utils.Triplet;
import view.generator.Stopper;
import view.generator.boxes.HBoxCriterion;
import view.generator.boxes.HBoxSolverCriterion;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import org.chocosolver.solver.constraints.Constraint;


public class GeneralModel {
    private Solver chocoSolver;
    private SolverResults solverResults;

    private final GeneratorRun generatorRun = new GeneratorRun(this);

    private final HashMap<String, Integer> variablesDegrees = new HashMap<>();

    // Dans la section des déclarations de champs de GeneralModel.java

    /**
     * Stores pairs of base hexagon indices that can undergo a C6-C6 -> C5-C7 transformation.
     */
    private ArrayList<Couple<Integer, Integer>> potentialTransformationSites;

    /**
     * Choco BoolVar array indicating if a transformation is active for each potential site.
     */
    private BoolVar[] transformationActiveVars;

    /**
     * Choco IntVar representing the number of C5 cycles to be generated via transformations.
     * This will be constrained by PentagonNumberConstraint.
     */
    private IntVar nbPentagonsVar;

    /**
     * Choco IntVar representing the number of C7 cycles to be generated via transformations.
     * This will be constrained by HeptagonNumberConstraint.
     */
    private IntVar nbHeptagonsVar;

    private final int nbMaxHexagons;

    private int[][] neighborIndices;

    //private final ArrayList<Couple<Integer, Integer>> outterHexagons = new ArrayList<>();
    private final ArrayList<Integer> outterHexagonsIndexes = new ArrayList<>();

    private ArrayList<ArrayList<Integer>> neighborGraphOutterHexagons;
    //private int indexOutterHexagon;

    private Couple<Integer, Integer>[] coordsCorrespondance;

    /*
     * Parameters
     */

    // Don't use for regular solving
    private final boolean applySymmetriesConstraints;

    private final int nbCrowns;
    //private int nbHexagons;
    private final int diameter;

    //private int nbEdges;
    private int nbClausesLexLead = 0;

    boolean verbose = false;

    /*
     * Constraint programming variables
     */

    private final Model chocoModel = new Model("Benzenoides");

    private Node[] nodesRefs;

    private int[][] hexagonIndicesMatrix;
    private int[] hexagonSparseIndicesTab;
    private int[] hexagonCompactIndicesTab;

    private int nbHexagonsCoronenoid;

    private int[][] sideSharing;
    private int[][] adjacencyMatrixWithOutterHexagons;

    /***
     * Choco solver variables
     */
    private UndirectedGraph GUB;

    private UndirectedGraphVar benzenoidGraphVar;
    private BoolVar[] hexBoolVars;
    private BoolVar[] benzenoidVerticesBVArray;
    private BoolVar[][] benzenoidEdges;

    private final ArrayList<Variable> variables = new ArrayList<>();

    private IntVar nbVertices;
    //private BoolVar[] edges;

    private BoolVar[] nbHexagonsReifies;
    //private IntVar graphDiameter;

    private final ArrayList<ArrayList<Integer>> nogoods = new ArrayList<>();

    private final SimpleIntegerProperty nbTotalSolutions = new SimpleIntegerProperty(0);
    private int indexSolution;


    /*
     * Properties
     */

    private final ModelPropertySet modelPropertySet;
    private static final SolverPropertySet solverPropertySet = new SolverPropertySet();


    private boolean isInTestMode = false;

    /*
     * Constructors
     */

    public GeneralModel(ModelPropertySet modelPropertySet) {
        this.modelPropertySet = modelPropertySet;
        nbMaxHexagons = modelPropertySet.computeHexagonNumberUpperBound();
        nbCrowns = modelPropertySet.computeNbCrowns();
        diameter = (2 * nbCrowns) - 1;
        applySymmetriesConstraints = modelPropertySet.symmetryConstraintsAppliable();
        initialize();
    }

    public GeneralModel(ModelPropertySet modelPropertySet, int nbCrowns) {
        this.modelPropertySet = modelPropertySet;
        nbMaxHexagons = modelPropertySet.computeHexagonNumberUpperBound();
        this.nbCrowns = nbCrowns;
        diameter = (2 * nbCrowns) - 1;
        applySymmetriesConstraints = modelPropertySet.symmetryConstraintsAppliable();
        initialize();
    }


    /*
     * Initialization methods
     */

    private void initialize() {
        nbHexagonsCoronenoid = 1 + (diameter + 1) * (diameter - 1) * 3 / 4;
        //hexagonIndices = initializeHexagonIndices(diameter);
        initializeVariables();
        initializeConstraints();
        initializeTransformationFramework(); // ADDED
        buildNodesRefs();
        System.out.print("");
    }


    private void initializeVariables() {
        System.out.println("00 " + nbMaxHexagons);
        nbHexagonsReifies = new BoolVar[nbMaxHexagons + 1];
        System.out.println("1");
        hexagonIndicesMatrix = buildHexagonIndices();
        hexagonSparseIndicesTab = buildHexagonSparseIndices(hexagonIndicesMatrix, diameter, nbHexagonsCoronenoid);
        hexagonCompactIndicesTab = buildHexagonCompactIndices(hexagonSparseIndicesTab, diameter);
        System.out.println("2");
        UndirectedGraph GLB = BoundsBuilder.buildGLB2(this);
        GUB = BoundsBuilder.buildGUB2(this);
        //indexOutterHexagon = diameter * diameter;
        System.out.println("3");
        sideSharing = buildAdjacencyMatrix();


        benzenoidGraphVar = chocoModel.graphVar("g", GLB, GUB);

        System.out.println("4");
        buildBenzenoidVertices();
        System.out.println("5");
        buildBenzenoidEdges();
        System.out.println("6");
        buildCoordsCorrespondance();
        System.out.println("7");
        buildNeighborIndices();

        nbVertices = chocoModel.intVar("nbVertices", 1, nbHexagonsCoronenoid);
        //graphDiameter = chocoModel.intVar("diameter", 0, diameter);

    }

    // Dans GeneralModel.java, à l'intérieur ou après initializeVariables()

    private void initializeTransformationFramework() {
        // Initialiser nbPentagonsVar et nbHeptagonsVar (sera fait par les classes de contraintes)
        // Par exemple, dans PentagonNumberConstraint.buildVariables():
        // generalModel.setNbPentagonsVar(chocoModel.intVar("nb_pentagons", 0, maxPotentialTransformations));
        // Idem pour nbHeptagonsVar.

        identifyPotentialTransformationSites(); // Vous devez implémenter cette méthode

        if (potentialTransformationSites != null && !potentialTransformationSites.isEmpty()) {
            transformationActiveVars = chocoModel.boolVarArray("transformationSite", potentialTransformationSites.size());
        } else {
            transformationActiveVars = new BoolVar[0]; // Aucun site potentiel
        }

        // S'assurer que nbPentagonsVar et nbHeptagonsVar sont bien initialisés par les contraintes correspondantes.
        // Si ce n'est pas le cas, initialisez-les ici avec une borne supérieure appropriée.
        // Par exemple, le nombre maximal de C5/C7 ne peut pas dépasser le nombre de sites.
        if (this.nbPentagonsVar == null) {
            int maxC5 = (potentialTransformationSites != null) ? potentialTransformationSites.size() : 0;
            this.nbPentagonsVar = chocoModel.intVar("nb_pentagons_default", 0, maxC5);
        }
        if (this.nbHeptagonsVar == null) {
            int maxC7 = (potentialTransformationSites != null) ? potentialTransformationSites.size() : 0;
            this.nbHeptagonsVar = chocoModel.intVar("nb_heptagons_default", 0, maxC7);
        }
    }

    /**
     * Identifies all pairs of adjacent base hexagons that can serve as potential
     * sites for a C6-C6 -> C5-C7 transformation.
     * This method needs to be implemented based on your grid structure.
     */
    private void identifyPotentialTransformationSites() {
        potentialTransformationSites = new ArrayList<>();
        // Parcourir la grille d'hexagones (hexagonIndicesMatrix)
        // Pour chaque hexagone h1, trouver ses voisins h2.
        // Si la paire (h1, h2) est un site valide pour la transformation, l'ajouter.
        // Attention à ne pas ajouter deux fois la même paire (ex: (h1,h2) et (h2,h1)).

        // Exemple de logique (à adapter à votre structure hexagonIndicesMatrix/sideSharing)
        if (hexagonIndicesMatrix == null || sideSharing == null) {
            System.err.println("Hexagon grid not initialized before identifying transformation sites.");
            return;
        }

        for (int i = 0; i < nbHexagonsCoronenoid; i++) { // Utilise le nombre d'hexagones dans la grille de base
            int h1_sparse = getHexagonSparseIndex(i); // Obtenir l'index "sparse" si nécessaire
            if (h1_sparse == -1) continue;

            for (int j = i + 1; j < nbHexagonsCoronenoid; j++) {
                int h2_sparse = getHexagonSparseIndex(j);
                if (h2_sparse == -1) continue;

                // Vérifier si h1_sparse et h2_sparse sont adjacents dans la grille de base
                // sideSharing utilise des index "sparse"
                if (sharesSide(h1_sparse, h2_sparse)) {
                    potentialTransformationSites.add(new Couple<>(h1_sparse, h2_sparse));
                }
            }
        }
        System.out.println("Identified " + potentialTransformationSites.size() + " potential transformation sites.");
    }

    // Méthode utilitaire pour obtenir l'index "sparse" à partir de l'index compact de hexBoolVars
// Ceci est un placeholder, adaptez-le à votre structure d'indexation.
    public int getHexagonSparseIndex(int compactIndex) {
        if (hexagonSparseIndicesTab != null && compactIndex < hexagonSparseIndicesTab.length) {
            return hexagonSparseIndicesTab[compactIndex];
        }
        // Gérer le cas où l'index n'est pas valide ou si la table n'est pas initialisée
        System.err.println("Warning: Could not get sparse index for compact index " + compactIndex);
        return -1; // Ou lancez une exception
    }


    //Getter pour les variables de transformation (nécessaire pour les contraintes)
    public BoolVar[] getTransformationActiveVars() {
        return transformationActiveVars;
    }

    public ArrayList<Couple<Integer, Integer>> getPotentialTransformationSites() {
        return potentialTransformationSites;
    }

    // Setters pour que les classes de contraintes puissent initialiser ces IntVars
    public void setNbPentagonsVar(IntVar nbPentagonsVar) {
        this.nbPentagonsVar = nbPentagonsVar;
    }

    public void setNbHeptagonsVar(IntVar nbHeptagonsVar) {
        this.nbHeptagonsVar = nbHeptagonsVar;
    }

    public IntVar getNbPentagonsVar() {
        return nbPentagonsVar;
    }

    public IntVar getNbHeptagonsVar() {
        return nbHeptagonsVar;
    }

    private int[] buildHexagonSparseIndices(int[][] hexagonIndices, int diameter, int nbHexagonsCoronenoid) {
        int [] hexagonSparseIndices = new int[nbHexagonsCoronenoid];
        int hexagonSparseIndex = 0;

        for (int lineIndex = 0; lineIndex < diameter; lineIndex++) {
            for (int columnIndex = 0; columnIndex < diameter; columnIndex++) {
                if (validHexagonIndex(lineIndex, columnIndex)) {
                    hexagonSparseIndices[hexagonSparseIndex] = hexagonIndices[lineIndex][columnIndex];
                    hexagonSparseIndex++;
                }
            }
        }
        return hexagonSparseIndices;
    }

    private int[] buildHexagonCompactIndices(int[] hexagonSparseIndices, int diameter) {
        int [] hexagonCompactIndices = new int[diameter * diameter];
        Arrays.fill(hexagonCompactIndices, -1);
        for (int i = 0; i < hexagonSparseIndices.length; i++)
            hexagonCompactIndices[hexagonSparseIndices[i]] = i;
        return hexagonCompactIndices;
    }

    private void initializeConstraints() {
        chocoModel.connected(benzenoidGraphVar).post();
        //chocoModel.diameter(benzenoidGraphVar,graphDiameter).post();
        //chocoModel.arithm(graphDiameter, "<", nbVertices).post();

        ConstraintBuilder.postFillNodesConnection(this);
        ConstraintBuilder.postNoHolesOfSize1Constraint(this);
        chocoModel.nbNodes(benzenoidGraphVar, nbVertices).post();
        if (applySymmetriesConstraints)
            nbClausesLexLead = ConstraintBuilder.postSymmetryBreakingConstraints(this);
    }

    public void addVariable(Variable variable) {
        variables.add(variable);
    }

    public void addVariable(Variable... variables) {
        for (Variable variable : variables)
            addVariable(variable);
    }

    /***
     * Apply the model property constraint to the model
     * @param modelProperty : the model property
     */
    public void applyModelConstraint(ModelProperty modelProperty) {
        modelProperty.getConstraint().build(this, modelProperty.getExpressions());
        postTransformationConstraints();
    }
    // Dans generator.GeneralModel.java

    private void postTransformationAtomExistsConstraints() {
        if (transformationActiveVars == null || potentialTransformationSites == null) { // Ajout d'une vérification
            System.err.println("Transformation framework non initialisé (transformationActiveVars ou potentialTransformationSites est null).");
            return;
        }

        for (int k = 0; k < potentialTransformationSites.size(); k++) {

            // Vérification de sécurité pour la taille de transformationActiveVars
            if (k >= transformationActiveVars.length) {
                System.err.println("Index k (" + k + ") hors limites pour transformationActiveVars (taille: " + transformationActiveVars.length + ")");
                continue;
            }

            Couple<Integer, Integer> site = potentialTransformationSites.get(k);
            int h1_sparse_idx = site.getX(); // Index "sparse"
            int h2_sparse_idx = site.getY(); // Index "sparse"

            int h1_compact_idx = getHexagonCompactIndex(h1_sparse_idx);
            int h2_compact_idx = getHexagonCompactIndex(h2_sparse_idx);

            if (h1_compact_idx == -1 || h2_compact_idx == -1) {
                // System.err.println("Index compact invalide pour le site (sparse): " + h1_sparse_idx + ", " + h2_sparse_idx);
                continue;
            }

            // Assurer que les index compacts sont valides pour hexBoolVars
            if (h1_compact_idx >= hexBoolVars.length || h2_compact_idx >= hexBoolVars.length) {
                System.err.println("Index compact hors limites pour hexBoolVars. h1_compact_idx: " + h1_compact_idx + ", h2_compact_idx: " + h2_compact_idx + ", taille hexBoolVars: " + hexBoolVars.length);
                continue;
            }

            // Définir la condition : la transformation au site k est active
            Constraint condition = getProblem().arithm(transformationActiveVars[k], "=", 1);

            // Définir la conséquence : les hexagones de base correspondants n'existent pas (en tant que C6)
            Constraint consequence = getProblem().and(
                    getProblem().arithm(hexBoolVars[h1_compact_idx], "=", 0),
                    getProblem().arithm(hexBoolVars[h2_compact_idx], "=", 0)
            );

            // Appliquer la contrainte ifThen. Choco la poste ou l'enregistre.
            getProblem().ifThen(condition, consequence);
            // PAS DE .post() ici car ifThen semble retourner void dans votre configuration.
        }
    }

    /***
     * Apply all the model constraints to the model
     */

    private void applyModelConstraints() {
        for (Property modelProperty : modelPropertySet)
            if (modelPropertySet.has(modelProperty.getId())) {
                applyModelConstraint((ModelProperty) modelProperty);
            }

        // Maintenant que nbPentagonsVar et nbHeptagonsVar devraient être initialisés par leurs contraintes respectives,
        // postez les contraintes de transformation.
        postTransformationConstraints(); // Assurez-vous que cette méthode existe et appelle les sous-méthodes
        if (!modelPropertySet.has("symmetry") && !modelPropertySet.has("rectangle") && !modelPropertySet.has("rhombus"))
            ConstraintBuilder.postBordersConstraints(this);
    }




    /*
     * Solving methods
     */

    public void displaySolution(Solver solver) {

        for (int index = 0; index < hexBoolVars.length; index++) {
            if (hexBoolVars[index].getValue() == 1)
                System.out.print(index + " ");
        }

        System.out.println();

        for (Variable x : variables)
            if ("XI".equals(x.getName()))
                System.out.println(x.getName() + " = " + (double) ((((IntVar) x).getValue())) / 100);
            else
                System.out.println(x.getName() + " = " + (((IntVar) x).getValue()));

        System.out.println(solver.getDecisionPath());

        if (!verbose)
            System.out.println();

        System.out.println(this.getProblem().getSolver().getDecisionPath());
        System.out.println(this.getProblem().getSolver().getFailCount() + " fails");
    }

    public String buildDescription(int index) {

        StringBuilder builder = new StringBuilder();
        builder.append("solution ").append(index).append("\n");

        for (Variable x : variables)

            if ("XI".equals(x.getName())) {

                double value = (double) ((((IntVar) x).getValue())) / 100;
                NumberFormat formatter = new DecimalFormat("#0.00");
                builder.append(x.getName()).append(" = ").append(formatter.format(value).replace(",", ".")).append("\n");

            } else
                builder.append(x.getName()).append(" = ").append(((IntVar) x).getValue()).append("\n");

        return builder.toString();
    }

    @SuppressWarnings("unchecked")
    private void buildCoordsCorrespondance() {

        coordsCorrespondance = new Couple[diameter * diameter];

        for (int lineIndex = 0; lineIndex < diameter; lineIndex++) {
            for (int columnIndex = 0; columnIndex < diameter; columnIndex++) {

                if (validHexagonIndex(lineIndex, columnIndex))
                    coordsCorrespondance[hexagonIndicesMatrix[lineIndex][columnIndex]] = new Couple<>(lineIndex, columnIndex);
            }
        }
    }

    private Pattern convertToPattern() {

        ArrayList<Integer> hexagonsSolutions = new ArrayList<>();

        int[] correspondance = new int[diameter * diameter];

        Arrays.fill(correspondance, -1);

        for (int index = 0; index < benzenoidVerticesBVArray.length; index++) {
            if (benzenoidVerticesBVArray[index] != null && benzenoidVerticesBVArray[index].getValue() == 1) {
                    hexagonsSolutions.add(index);
                    correspondance[index] = hexagonsSolutions.size() - 1;
            }
        }

        int nbNodes = hexagonsSolutions.size();

        /*
         * nodes
         */

        Node[] nodes = new Node[nbNodes];

        for (int i = 0; i < hexagonsSolutions.size(); i++) {

            int hexagon = hexagonsSolutions.get(i);

            Couple<Integer, Integer> couple = coordsCorrespondance[hexagon];

            nodes[i] = new Node(couple.getY(), couple.getX(), i);
        }

        /*
         * matrix
         */

        int[][] matrix = new int[nbNodes][nbNodes];
        int[][] neighbors = new int[nbNodes][6];

        for (int nodeIndex = 0; nodeIndex < nbNodes; nodeIndex++)
            Arrays.fill(neighbors[nodeIndex], -1);

        for (int nodeIndex = 0; nodeIndex < nbNodes; nodeIndex++) {

            int u = hexagonsSolutions.get(nodeIndex);
            Node n1 = nodes[nodeIndex];

            for (int j = (nodeIndex + 1); j < nbNodes; j++) {

                int v = hexagonsSolutions.get(j);
                Node n2 = nodes[j];

                if (sideSharing[u][v] == 1) {

                    // Setting matrix
                    matrix[nodeIndex][j] = 1;
                    matrix[j][nodeIndex] = 1;

                    // Setting neighbors
                    int x1 = n1.getX();
                    int y1 = n1.getY();
                    int x2 = n2.getX();
                    int y2 = n2.getY();

                    if (x2 == x1 && y2 == y1 - 1) {
                        neighbors[correspondance[u]][0] = correspondance[v];
                        neighbors[correspondance[v]][3] = correspondance[u];
                    } else if (x2 == x1 + 1 && y2 == y1) {
                        neighbors[correspondance[u]][1] = correspondance[v];
                        neighbors[correspondance[v]][4] = correspondance[u];
                    } else if (x2 == x1 + 1 && y2 == y1 + 1) {
                        neighbors[correspondance[u]][2] = correspondance[v];
                        neighbors[correspondance[v]][5] = correspondance[u];
                    } else if (x2 == x1 && y2 == y1 + 1) {
                        neighbors[correspondance[u]][3] = correspondance[v];
                        neighbors[correspondance[v]][0] = correspondance[u];
                    } else if (x2 == x1 - 1 && y2 == y1) {
                        neighbors[correspondance[u]][4] = correspondance[v];
                        neighbors[correspondance[v]][1] = correspondance[u];
                    } else if (x2 == x1 - 1 && y2 == y1 - 1) {
                        neighbors[correspondance[u]][5] = correspondance[v];
                        neighbors[correspondance[v]][2] = correspondance[u];
                    }
                }
            }
        }

        /*
         * Label
         */

        PatternLabel[] labels = new PatternLabel[nbNodes];
        Arrays.fill(labels, PatternLabel.POSITIVE);
        return new Pattern(matrix, labels, nodes, null, neighbors, 0);
    }

    private void recordNoGoods() {

        ArrayList<Integer> vertices = new ArrayList<>();
        for (int i = 0; i < hexBoolVars.length; i++)
            if (hexBoolVars[i].getValue() == 1)
                vertices.add(i);

        int center = hexagonCompactIndicesTab[hexagonIndicesMatrix[(diameter - 1) / 2][(diameter - 1) / 2]];

        Solution solution = new Solution(nodesRefs, hexagonCompactIndicesTab, hexagonSparseIndicesTab, hexagonIndicesMatrix,
                center, nbCrowns, vertices);

        NoGoodRecorder noGoodRecorder;// = null;

        if (!modelPropertySet.has("symmetry")) {
            ArrayList<Integer> topBorder = new ArrayList<>();
            for (int columnIndex = 0; columnIndex < diameter; columnIndex++) {
                if (validHexagonIndex(0, columnIndex))
                    topBorder.add(hexagonCompactIndicesTab[hexagonIndicesMatrix[0][columnIndex]]);
            }
            ArrayList<Integer> leftBorder = new ArrayList<>();
            for (int lineIndex = 0; lineIndex < diameter; lineIndex++) {
                for (int columnIndex = 0; columnIndex < diameter; columnIndex++) {
                    if (validHexagonIndex(lineIndex, columnIndex)) {
                        leftBorder.add(hexagonCompactIndicesTab[hexagonIndicesMatrix[lineIndex][columnIndex]]);
                        break;
                    }
                }
            }
            solution.setPattern(convertToPattern());
            noGoodRecorder = new NoGoodBorderRecorder(this, solution, topBorder, leftBorder);
        } else {

            noGoodRecorder = new NoGoodNoneRecorder(this, solution);
            //TODO changer tests
            if (Objects.equals(((ParameterizedExpression) modelPropertySet.getById("symmetry").getExpressions().get(0)).getOperator(), "SYMM_MIRROR"))
                noGoodRecorder = new NoGoodHorizontalAxisRecorder(this, solution);
            else if (Objects.equals(((ParameterizedExpression) modelPropertySet.getById("symmetry").getExpressions().get(0)).getOperator(), "SYMM_VERTICAL"))
                //noGoodRecorder = new NoGoodHorizontalAxisRecorder(this, solution);
                noGoodRecorder = new NoGoodVerticalAxisRecorder(this, solution);

            else {
                if (modelPropertySet.has("pattern"))
                    noGoodRecorder = new NoGoodUniqueRecorder(this, solution);
            }
        }
        solution.setPattern(convertToPattern());
        noGoodRecorder = new NoGoodAllRecorder(this, solution);
        noGoodRecorder.record();

    }

    // Dans generator.GeneralModel.java

    public SolverResults solve() {
        // Appliquer les contraintes de base du modèle et celles définies par l'utilisateur
        applyModelConstraints(); // S'assure que nbPentagonsVar et nbHeptagonsVar sont initialisés
        // et que les contraintes de transformation sont postées.

        chocoSolver = chocoModel.getSolver(); // Obtenir le solveur après que toutes les variables et contraintes sont postées.

        // Définir la stratégie de recherche
        // Inclure hexBoolVars et transformationActiveVars si ce dernier est non nul et non vide
        IntVar[] allDecisionVars;
        if (transformationActiveVars != null && transformationActiveVars.length > 0) {
            allDecisionVars = new IntVar[hexBoolVars.length + transformationActiveVars.length];
            System.arraycopy(hexBoolVars, 0, allDecisionVars, 0, hexBoolVars.length);
            System.arraycopy(transformationActiveVars, 0, allDecisionVars, hexBoolVars.length, transformationActiveVars.length);
        } else {
            allDecisionVars = hexBoolVars;
        }
        // Appliquer une stratégie de recherche par défaut. Elle peut être modifiée par des contraintes spécifiques.
        chocoSolver.setSearch(new IntStrategy(allDecisionVars, new FirstFail(chocoModel), new IntDomainMin()));


        // Permettre aux contraintes individuelles de modifier la stratégie de recherche si nécessaire.
        // Attention : cela pourrait surcharger la stratégie définie ci-dessus.
        // Il faut s'assurer que cela reste cohérent.
        for (Property modelProperty : modelPropertySet) {
            if (modelProperty.hasExpressions()) {
                ((ModelProperty) modelProperty).getConstraint().changeSolvingStrategy();
            }
        }

        // Appliquer les propriétés du solveur (limites de temps, nombre de solutions)
        for (Property solverProperty : solverPropertySet) {
            if (solverProperty.hasExpressions()) {
                ((SolverProperty) solverProperty).getSpecifier().apply(chocoSolver, solverProperty.getExpressions().get(0));
            }
        }

        solverResults = new SolverResults();
        indexSolution = 0;

        long begin = System.currentTimeMillis();

        chocoSolver.limitSearch(() -> Stopper.STOP); // Permet d'arrêter la recherche depuis l'extérieur
        Stopper.STOP = false; // Réinitialiser le stopper

        while (chocoSolver.solve() && !generatorRun.isPaused()) {

            // À ce stade, Choco a trouvé une affectation pour hexBoolVars et transformationActiveVars
            // qui satisfait toutes les contraintes.

            ArrayList<Integer> verticesSolution = buildVerticesSolution(); // Basé sur hexBoolVars
            String description = buildDescription(indexSolution);

            // LA CONSTRUCTION DE LA MOLECULE DEVIENT PLUS COMPLEXE ICI
            // Vous devez maintenant interpréter `verticesSolution` ET `transformationActiveVars`
            // pour construire la géométrie atomique réelle avec les cycles C5, C6, C7.
            // L'objet `Benzenoid` retourné par `buildMolecule` doit refléter cette nouvelle structure.
            // Pour l'instant, on continue avec l'ancienne méthode, mais il faudra la réviser.

            // ----- DÉBUT SECTION À REVOIR PROFONDÉMENT -----
            Benzenoid molecule = Benzenoid.buildMolecule(description, nbCrowns, indexSolution, verticesSolution);
            // Après la création de la molécule de base (hexagonale), il faudrait :
            // 1. Identifier les transformations C5/C7 actives (via transformationActiveVars).
            // 2. Modifier la structure de 'molecule' pour refléter ces transformations
            //    (supprimer/ajouter des atomes, changer la connectivité).
            // 3. Recalculer les faces (cycles) de la molécule modifiée (C5, C6, C7).
            // 4. Mettre à jour molecule.setExplicitCycleTypes(...)
            // ----- FIN SECTION À REVOIR PROFONDÉMENT -----

            if (molecule.respectPostProcessing(modelPropertySet)) { // Ce filtre doit aussi être conscient des C5/C7
                solverResults.addMolecule(molecule);
                if (inTestMode()) {
                    nbTotalSolutions.set(nbTotalSolutions.get() + 1);
                } else {
                    Platform.runLater(() -> nbTotalSolutions.set(nbTotalSolutions.get() + 1));
                }

                recordNoGoods(); // Cette méthode doit aussi potentiellement être adaptée

                // BenzenoidSolution est probablement pour une représentation purement hexagonale.
                // À revoir si la structure interne change.
                BenzenoidSolution solverSolution = new BenzenoidSolution(GUB, nbCrowns,
                        chocoModel.getName() + indexSolution, hexagonSparseIndicesTab);

                solverResults.addSolution(solverSolution, description, nbCrowns);
                solverResults.addVerticesSolution(verticesSolution); // Conserver la solution Choco brute pour le débogage

                displaySolution(chocoSolver); // Pour le débogage

                if (verbose) {
                    System.out.println("NO-GOOD");
                    for (ArrayList<Integer> ng : nogoods) {
                        for (Integer v : ng)
                            System.out.print(v + " ");
                        System.out.println();
                    }
                }
                indexSolution++;
            }
        }

        long end = System.currentTimeMillis();
        long time = end - begin;

        solverResults.setTime(time);
        solverResults.setNbTotalSolution(nbTotalSolutions.get());
        solverResults.setSolver(chocoSolver);

        System.out.println(nbCrowns + " crowns");
        System.out.println(nogoods.size() + " no-good clauses");
        System.out.println(nbClausesLexLead + " lex-lead clauses");
        chocoSolver.printStatistics();

        solverResults.setNogoodsFragments();
        System.out.println("------");
        displayDegrees();
        return solverResults;
    }

    private boolean inTestMode() {
        return isInTestMode;
    }


    private ArrayList<Integer> buildVerticesSolution() {
        ArrayList<Integer> verticesSolution = new ArrayList<>();

        for (BoolVar benzenoidVertex : benzenoidVerticesBVArray) {

            if (benzenoidVertex != null) {
                verticesSolution.add(benzenoidVertex.getValue());
            } else
                verticesSolution.add(0);
        }
        return verticesSolution;
    }

    private int[][] buildAdjacencyMatrix() {

        int nbEdges = 0;
        int [][] adjacencyMatrix = new int[diameter * diameter][diameter * diameter];

        for (int lineIndex = 0; lineIndex < diameter; lineIndex++) {
            for (int columnIndex = 0; columnIndex < diameter; columnIndex++) {

                if (validHexagonIndex(lineIndex, columnIndex)) {

                    int u = hexagonIndicesMatrix[lineIndex][columnIndex];

                    if (lineIndex > 0 && columnIndex > 0) {

                        int v = hexagonIndicesMatrix[lineIndex - 1][columnIndex - 1];
                        if (v != -1) {
                            if (adjacencyMatrix[u][v] == 0) {
                                adjacencyMatrix[u][v] = 1;
                                adjacencyMatrix[v][u] = 1;
                                nbEdges++;
                            }
                        }
                    }

                    if (columnIndex > 0) {

                        int v = hexagonIndicesMatrix[lineIndex][columnIndex - 1];
                        if (v != -1) {
                            if (adjacencyMatrix[u][v] == 0) {
                                adjacencyMatrix[u][v] = 1;
                                adjacencyMatrix[v][u] = 1;
                                nbEdges++;
                            }
                        }
                    }

                    if (lineIndex + 1 < diameter) {

                        int v = hexagonIndicesMatrix[lineIndex + 1][columnIndex];
                        if (v != -1) {
                            if (adjacencyMatrix[u][v] == 0) {
                                adjacencyMatrix[u][v] = 1;
                                adjacencyMatrix[v][u] = 1;
                                nbEdges++;
                            }
                        }
                    }

                    if (lineIndex + 1 < diameter && columnIndex + 1 < diameter) {

                        int v = hexagonIndicesMatrix[lineIndex + 1][columnIndex + 1];
                        if (v != -1) {
                            if (adjacencyMatrix[u][v] == 0) {
                                adjacencyMatrix[u][v] = 1;
                                adjacencyMatrix[v][u] = 1;
                                nbEdges++;
                            }
                        }
                    }

                    if (columnIndex + 1 < diameter) {

                        int v = hexagonIndicesMatrix[lineIndex][columnIndex + 1];
                        if (v != -1) {
                            if (adjacencyMatrix[u][v] == 0) {
                                adjacencyMatrix[u][v] = 1;
                                adjacencyMatrix[v][u] = 1;
                                nbEdges++;
                            }
                        }
                    }

                    if (lineIndex > 0) {

                        int v = hexagonIndicesMatrix[lineIndex - 1][columnIndex];
                        if (v != -1) {
                            if (adjacencyMatrix[u][v] == 0) {
                                adjacencyMatrix[u][v] = 1;
                                adjacencyMatrix[v][u] = 1;
                                nbEdges++;
                            }
                        }
                    }
                }
            }
        }
        return adjacencyMatrix;
    }

    /*public int getNbEdges() {
        return nbEdges;
    }*/

    public ArrayList<Integer> getOutterHexagonsIndexes() {
        return outterHexagonsIndexes;
    }

    public void buildNeighborGraphWithOutterHexagons(int order) {

        if (order == 0) {
            // buildNeighborIndices();
            neighborGraphOutterHexagons = new ArrayList<>();

            for (int[] ints : neighborIndices) {
                ArrayList<Integer> neighbors = new ArrayList<>();
                for (int neighborIndex = 0; neighborIndex < 6; neighborIndex++)
                    neighbors.add(ints[neighborIndex]);
                neighborGraphOutterHexagons.add(neighbors);
            }
        } else {

            ArrayList<Integer> hexagons = new ArrayList<>();
            ArrayList<Couple<Integer, Integer>> coords = new ArrayList<>();

            for (int lineIndex = 0; lineIndex < diameter; lineIndex++) {
                if (lineIndex == 0 || lineIndex == diameter - 1) {
                    for (int columnIndex = 0; columnIndex < diameter; columnIndex++) {
                        if (validHexagonIndex(lineIndex, columnIndex)) {
                            hexagons.add(hexagonIndicesMatrix[lineIndex][columnIndex]);
                            coords.add(new Couple<>(columnIndex, lineIndex));
                        }
                    }
                } else {

                    ArrayList<Integer> c1 = getFirstColumns(lineIndex, 1);
                    ArrayList<Integer> c2 = getLastColumns(lineIndex, 1);

                    for (int i = 0; i < c1.size(); i++) {

                        hexagons.add(hexagonIndicesMatrix[lineIndex][c1.get(i)]);
                        coords.add(new Couple<>(c1.get(i), lineIndex));
                        hexagons.add(hexagonIndicesMatrix[lineIndex][c2.get(i)]);
                        coords.add(new Couple<>(c2.get(i), lineIndex));
                    }
                }
            }

            // buildNeighborIndices();

            neighborGraphOutterHexagons = new ArrayList<>();

            for (int[] ints : neighborIndices) {
                ArrayList<Integer> neighbors = new ArrayList<>();
                for (int j = 0; j < 6; j++)
                    neighbors.add(ints[j]);
                neighborGraphOutterHexagons.add(neighbors);
            }

            ArrayList<Triplet<Integer, Integer, Integer>> coordsOutterHexagons = new ArrayList<>();

            for (int i = 0; i < hexagons.size(); i++) {

                int hexagon = hexagons.get(i);
                Couple<Integer, Integer> coord = coords.get(i);

                int[] neighbors = neighborIndices[hexagon];

                for (HexNeighborhood neighbor : HexNeighborhood.values()) {
                    if (neighbors[neighbor.getIndex()] == -1) {
                        int x = coord.getX() + neighbor.dx();
                        int y = coord.getY() + neighbor.dy();

                        int indexOutter = -1;

                        for (Triplet<Integer, Integer, Integer> coordOutter : coordsOutterHexagons)
                            if (coordOutter.getX() == x && coordOutter.getY() == y)
                                indexOutter = coordOutter.getZ();

                        if (indexOutter == -1) {

                            indexOutter = neighborGraphOutterHexagons.size();

                            coordsOutterHexagons.add(new Triplet<>(x, y, indexOutter));

                            ArrayList<Integer> newNeighbor = new ArrayList<>();
                            for (int k = 0; k < 6; k++)
                                newNeighbor.add(-1);

                            neighborGraphOutterHexagons.add(newNeighbor);
                        }

                        neighborGraphOutterHexagons.get(hexagon).set(neighbor.getIndex(), indexOutter);
                        neighborGraphOutterHexagons.get(indexOutter).set((neighbor.getIndex() + 3) % 6, hexagon);
                    }
                }
            }

            for (int i = 0; i < coordsOutterHexagons.size(); i++) {

                Triplet<Integer, Integer, Integer> coord1 = coordsOutterHexagons.get(i);

                int x1 = coord1.getX();
                int y1 = coord1.getY();

                int index1 = coord1.getZ();

                for (int j = 0; j < coordsOutterHexagons.size(); j++) {

                    if (i != j) {

                        Triplet<Integer, Integer, Integer> coord2 = coordsOutterHexagons.get(j);

                        int x2 = coord2.getX();
                        int y2 = coord2.getY();

                        int index2 = coord2.getZ();

                        if (x2 == x1 && y2 == y1 - 1) {
                            neighborGraphOutterHexagons.get(index1).set(0, index2);
                            neighborGraphOutterHexagons.get(index2).set(3, index1);
                        } else if (x2 == x1 + 1 && y2 == y1) {
                            neighborGraphOutterHexagons.get(index1).set(1, index2);
                            neighborGraphOutterHexagons.get(index2).set(4, index1);
                        } else if (x2 == x1 + 1 && y2 == y1 + 1) {
                            neighborGraphOutterHexagons.get(index1).set(2, index2);
                            neighborGraphOutterHexagons.get(index2).set(5, index1);
                        } else if (x2 == x1 && y2 == y1 + 1) {
                            neighborGraphOutterHexagons.get(index1).set(3, index2);
                            neighborGraphOutterHexagons.get(index2).set(0, index1);
                        } else if (x2 == x1 - 1 && y2 == y1) {
                            neighborGraphOutterHexagons.get(index1).set(4, index2);
                            neighborGraphOutterHexagons.get(index2).set(1, index1);
                        } else if (x2 == x1 - 1 && y2 == y1 - 1) {
                            neighborGraphOutterHexagons.get(index1).set(5, index2);
                            neighborGraphOutterHexagons.get(index2).set(2, index1);
                        }
                    }

                }
            }
        }
    }

    private ArrayList<Integer> getFirstColumns(int lineIndex, int order) {

        ArrayList<Integer> columns = new ArrayList<>();

        int nbColumns = 0;
        int columnIndex = 0;

        while (nbColumns < order && columnIndex < diameter) {

            if (validHexagonIndex(lineIndex, columnIndex)) {
                columns.add(columnIndex);
                nbColumns++;
            }

            columnIndex++;
        }

        return columns;
    }

    private ArrayList<Integer> getLastColumns(int lineIndex, int order) {

        ArrayList<Integer> columns = new ArrayList<>();

        int nbColumns = 0;
        int columnIndex = diameter - 1;

        while (nbColumns < order && columnIndex >= 0) {

            if (validHexagonIndex(lineIndex, columnIndex)) {
                columns.add(columnIndex);
                nbColumns++;
            }

            columnIndex--;
        }

        return columns;
    }

    private void buildNeighborIndices() {

        neighborIndices = new int[diameter * diameter][6];

        for (int[] ints : neighborIndices) {
            Arrays.fill(ints, -1);
        }

        for (int lineIndex = 0; lineIndex < hexagonIndicesMatrix.length; lineIndex++) {
            for (int columnIndex = 0; columnIndex < hexagonIndicesMatrix[lineIndex].length; columnIndex++) {

                if (validHexagonIndex(lineIndex, columnIndex)) {

                    int hexagonIndex = hexagonIndicesMatrix[lineIndex][columnIndex];

                    // Top-Right
                    if (lineIndex > 0)
                        neighborIndices[hexagonIndex][0] = hexagonIndicesMatrix[lineIndex - 1][columnIndex];

                    // Right
                    if (columnIndex < hexagonIndicesMatrix[lineIndex].length - 1)
                        neighborIndices[hexagonIndex][1] = hexagonIndicesMatrix[lineIndex][columnIndex + 1];

                    // Bottom-Right
                    if (lineIndex < hexagonIndicesMatrix[lineIndex].length - 1 && columnIndex < hexagonIndicesMatrix[lineIndex].length - 1)
                        neighborIndices[hexagonIndex][2] = hexagonIndicesMatrix[lineIndex + 1][columnIndex + 1];

                    // Bottom-Left
                    if (lineIndex < hexagonIndicesMatrix[lineIndex].length - 1)
                        neighborIndices[hexagonIndex][3] = hexagonIndicesMatrix[lineIndex + 1][columnIndex];

                    // Left
                    if (columnIndex > 0)
                        neighborIndices[hexagonIndex][4] = hexagonIndicesMatrix[lineIndex][columnIndex - 1];

                    // Top-Left
                    if (lineIndex > 0 && columnIndex > 0)
                        neighborIndices[hexagonIndex][5] = hexagonIndicesMatrix[lineIndex - 1][columnIndex - 1];
                }
            }
        }
    }

    public void buildAdjacencyMatrixWithOutterHexagons() {

        adjacencyMatrixWithOutterHexagons = new int[neighborGraphOutterHexagons.size()][neighborGraphOutterHexagons
                .size()];

        for (int nodeIndex = 0; nodeIndex < neighborGraphOutterHexagons.size(); nodeIndex++) {

            for (int neighborIndex = 0; neighborIndex < 6; neighborIndex++) {

                int v = neighborGraphOutterHexagons.get(nodeIndex).get(neighborIndex);

                if (v != -1) {
                    adjacencyMatrixWithOutterHexagons[nodeIndex][v] = 1;
                    adjacencyMatrixWithOutterHexagons[v][nodeIndex] = 1;
                }
            }
        }
    }

    public ArrayList<ArrayList<Integer>> getNeighborGraphOutterHexagons() {
        return neighborGraphOutterHexagons;
    }

    public int[][] getAdjacencyMatrixOutterHexagons() {
        return adjacencyMatrixWithOutterHexagons;
    }

    private void buildBenzenoidVertices() {

        hexBoolVars = new BoolVar[nbHexagonsCoronenoid];
        for (int i = 0; i < hexBoolVars.length; i++)
            hexBoolVars[i] = chocoModel.boolVar("vertex[" + i + "]");

        chocoModel.nodesChanneling(benzenoidGraphVar, hexBoolVars).post();

        benzenoidVerticesBVArray = new BoolVar[diameter * diameter];

        int index = 0;

        for (int lineIndex = 0; lineIndex < diameter; lineIndex++) {
            for (int columnIndex = 0; columnIndex < diameter; columnIndex++) {
                if (validHexagonIndex(lineIndex, columnIndex)) {

                    BoolVar hexBV = hexBoolVars[hexagonCompactIndicesTab[hexagonIndicesMatrix[lineIndex][columnIndex]]];
                    benzenoidVerticesBVArray[index] = hexBV;
                }
                index++;
            }
        }
    }

    private void buildBenzenoidEdges() {
        benzenoidEdges = new BoolVar[nbHexagonsCoronenoid][nbHexagonsCoronenoid];
        for (int lineIndex = 0; lineIndex < sideSharing.length; lineIndex++) {
            for (int columnIndex = (lineIndex + 1); columnIndex < sideSharing.length; columnIndex++) {
                if (sideSharing[lineIndex][columnIndex] == 1) {
                    int hexCompactIndex1 = hexagonCompactIndicesTab[lineIndex];
                    int hexCompactIndex2 = hexagonCompactIndicesTab[columnIndex];
                    BoolVar benzenoidEdgeBoolVar = chocoModel.boolVar("e_" + hexCompactIndex1 + "_" + hexCompactIndex2);
                    benzenoidEdges[hexCompactIndex1][hexCompactIndex2] = benzenoidEdgeBoolVar;
                    benzenoidEdges[hexCompactIndex2][hexCompactIndex1] = benzenoidEdgeBoolVar;
                    chocoModel.edgeChanneling(benzenoidGraphVar, benzenoidEdgeBoolVar, hexCompactIndex1, hexCompactIndex2).post();
                }
            }
        }
    }

    private int[][] buildHexagonIndices() {

        int [][] hexagonIndices = initializeHexagonIndices(diameter);

        int hexagonSparseIndex = 0;
        int centerIndex = (diameter - 1) / 2;

        int shift = diameter - nbCrowns;

        for (int lineIndex = 0; lineIndex < centerIndex; lineIndex++) {

            for (int columnIndex = 0; columnIndex < diameter - shift; columnIndex++) {
                hexagonIndices[lineIndex][columnIndex] = hexagonSparseIndex;
                hexagonSparseIndex++;
            }
            hexagonSparseIndex += shift;
            shift--;
        }

        for (int columnIndex = 0; columnIndex < diameter; columnIndex++) {
            hexagonIndices[centerIndex][columnIndex] = hexagonSparseIndex;
            hexagonSparseIndex++;
        }

        shift = 1;

        for (int lineIndex = centerIndex + 1; lineIndex < diameter; lineIndex++) {
            hexagonSparseIndex += shift;
            for (int columnIndex = shift; columnIndex < diameter; columnIndex++) {
                hexagonIndices[lineIndex][columnIndex] = hexagonSparseIndex;
                hexagonSparseIndex++;
            }
            shift++;
        }

        return hexagonIndices;
    }
    private int[][] initializeHexagonIndices(int diameter) {
        int[][] hexagonIndices = new int[diameter][diameter];
        for (int lineIndex = 0; lineIndex < diameter; lineIndex++) {
            Arrays.fill(hexagonIndices[lineIndex], -1);
        }
        return hexagonIndices;
    }
    public boolean validHexagonIndex(int lineIndex, int columnIndex) {
        return hexagonIndicesMatrix[lineIndex][columnIndex] != -1;
    }

    public int getNbHexagonsCoronenoid() {
        return nbHexagonsCoronenoid;
    }

    public int[] getHexagonCompactIndicesTab() {
        return hexagonCompactIndicesTab;
    }

    public int getHexagonCompactIndex(int i) {
        return hexagonCompactIndicesTab[i];
    }
    public BoolVar[] getNeighbors(int lineIndex, int columnIndex) {

        if (validHexagonIndex(lineIndex, columnIndex)) {

            BoolVar[] N = new BoolVar[6];

            for (int k = 0; k < 6; k++)
                N[k] = null;

            for(HexNeighborhood neighbor : HexNeighborhood.values()){
                int lineIndex2 = lineIndex + neighbor.dy();
                int columnIndex2 = columnIndex + neighbor.dx();
                if(lineIndex2 >= 0 && lineIndex2 <= diameter - 1 && columnIndex2 >= 0 && columnIndex2 <= diameter - 1)
                    if(validHexagonIndex(lineIndex2, columnIndex2))
                        N[neighbor.getIndex()] = benzenoidVerticesBVArray[hexagonIndicesMatrix[lineIndex2][columnIndex2]];
            }
            return N;
        }
        return null;
    }

    public BoolVar[] getHexBoolVars() {
        return hexBoolVars;
    }

    public BoolVar getHexBoolVar(int i) {
        return hexBoolVars[i];
    }

    public int[][] getNeighborIndices() {
        return neighborIndices;
    }

    public int getNbMaxHexagons() {
        return nbMaxHexagons;
    }

    public IntVar getNbVerticesVar() {
        return nbVertices;
    }

    public SolverResults getResultSolver() {
        return solverResults;
    }

    public GeneratorRun getGeneratorRun() {
        return generatorRun;
    }

    public void stop() {
        generatorRun.stop();
    }

    public boolean isPaused() {
        return generatorRun.isPaused();
    }

    public void pause() {
        generatorRun.pause();
    }

    public void resume() {

        while (chocoSolver.solve() && !generatorRun.isPaused()) {

            Platform.runLater(() -> nbTotalSolutions.set(nbTotalSolutions.get() + 1));

            recordNoGoods();

            BenzenoidSolution solution = new BenzenoidSolution(GUB, nbCrowns, chocoModel.getName() + indexSolution,
                    hexagonSparseIndicesTab);
            String description = buildDescription(indexSolution);
            solverResults.addSolution(solution, description, nbCrowns);

            ArrayList<BoolVar> presentHexagons = new ArrayList<>();
            ArrayList<Integer> verticesSolution = new ArrayList<>();

            for (BoolVar benzenoidVertex : benzenoidVerticesBVArray) {

                if (benzenoidVertex != null) {
                    verticesSolution.add(benzenoidVertex.getValue());

                    if (benzenoidVertex.getValue() == 1) {
                        presentHexagons.add(benzenoidVertex);
                    }

                } else
                    verticesSolution.add(0);
            }

            solverResults.addVerticesSolution(verticesSolution);

            displaySolution(chocoSolver);
            //System.out.println(solver.getDecisionPath());

            if (verbose) {

                System.out.println("NO-GOOD");

                for (ArrayList<Integer> ng : nogoods) {
                    for (Integer v : ng)
                        System.out.println(v + " ");
                }

                System.out.println();
            }

            indexSolution++;
        }

    }

    private void buildNodesRefs() {
        nodesRefs = new Node[nbHexagonsCoronenoid];

        for (int lineIndex = 0; lineIndex < diameter; lineIndex++) {
            for (int columnIndex = 0; columnIndex < diameter; columnIndex++) {
                if (validHexagonIndex(lineIndex, columnIndex)) {
                    int index = hexagonCompactIndicesTab[hexagonIndicesMatrix[lineIndex][columnIndex]];
                    nodesRefs[index] = new Node(lineIndex, columnIndex, index);
                }
            }
        }
    }

    private boolean isValid(Couple<Integer, Integer> coord, Pattern fragment, int index) {
        if (coord.getX() < 0 || coord.getX() >= diameter || coord.getY() < 0 || coord.getY() >= diameter
                || hexagonIndicesMatrix[coord.getY()][coord.getX()] == -1) {
            return fragment.getLabel(index) != PatternLabel.POSITIVE;
        }
        return true;
    }

    private int findIndex(ArrayList<Couple<Integer, Integer>> coords, Couple<Integer, Integer> coord) {
        for (int i = 0; i < coords.size(); i++)
            if (coords.get(i).equals(coord))
                return i;
        return -1;
    }

    @SuppressWarnings("unchecked")
    public PatternOccurences computeTranslations(Pattern pattern) {

        PatternOccurences fragmentOccurences = new PatternOccurences();

        /*
         * Trouver l'hexagone pr�sent du fragment le plus en haut � gauche
         */

        int minY = Integer.MAX_VALUE;
        for (Node node : pattern.getNodesRefs())
            if (node.getY() < minY)
                minY = node.getY();

        while (true) {
            boolean containsPresentHexagon = false;
            for (int i = 0; i < pattern.getNbNodes(); i++) {
                Node node = pattern.getNodesRefs()[i];
                if (node.getY() == minY && pattern.getLabel(i) == PatternLabel.POSITIVE)
                    containsPresentHexagon = true;
            }
            if (containsPresentHexagon)
                break;
            minY++;
        }

        int nodeIndex = -1;
        int minX = Integer.MAX_VALUE;
        for (int i = 0; i < pattern.getNbNodes(); i++) {
            Node node = pattern.getNodesRefs()[i];
            if (node.getY() == minY && node.getX() < minX && pattern.getLabel(i) == PatternLabel.POSITIVE) {
                minX = node.getX();
                nodeIndex = i;
            }
        }

        /*
         * Trouver les positions ou le fragment peut �tre plac�
         */

        for (int lineIndex = 0; lineIndex < diameter; lineIndex++) {
            for (int columnIndex = 0; columnIndex < diameter; columnIndex++) {
                int hexagonIndex = hexagonIndicesMatrix[lineIndex][columnIndex];
                if (hexagonIndex != -1) {
                    /*
                     * On place le fragment dans le coron�no�de de telle sorte que firstNode
                     * corresponde � hexagon
                     */

                    int[] checkedHexagons = new int[pattern.getNbNodes()];
                    Couple<Integer, Integer>[] coords = new Couple[pattern.getNbNodes()];

                    int candidat = nodeIndex;
                    checkedHexagons[nodeIndex] = 1;
                    coords[nodeIndex] = new Couple<>(columnIndex, lineIndex);

                    ArrayList<Integer> candidats = new ArrayList<>();

                    for (HexNeighborhood neighbor : HexNeighborhood.values()) {
                        if (pattern.getNeighbor(candidat, neighbor.getIndex()) != -1) {
                            int neighborIndex = pattern.getNeighbor(candidat, neighbor.getIndex());
                            candidats.add(neighborIndex);
                            coords[neighborIndex] = new Couple<>(columnIndex + neighbor.dx(), lineIndex + neighbor.dy());
                            checkedHexagons[neighborIndex] = 1;
                        }
                    }

                    while (candidats.size() > 0) {

                        candidat = candidats.get(0);

                        for (HexNeighborhood neighbor : HexNeighborhood.values()) {
                            if (pattern.getNeighbor(candidat, neighbor.getIndex()) != -1) {

                                int neighborIndex = pattern.getNeighbor(candidat, neighbor.getIndex());

                                if (checkedHexagons[neighborIndex] == 0) {
                                    candidats.add(neighborIndex);
                                    coords[neighborIndex] = new Couple<>(coords[candidat].getX() + neighbor.dx(), coords[candidat].getY() + neighbor.dy());
                                    checkedHexagons[neighborIndex] = 1;
                                }
                            }
                        }

                        candidats.remove(candidats.get(0));
                    }

                    /*
                     * On teste si le fragment obtenu est valide
                     */

                    boolean valid = true;
                    for (int i = 0; i < coords.length; i++) {
                        Couple<Integer, Integer> coord = coords[i];
                        if (!isValid(coord, pattern, i))
                            valid = false;
                    }
                    if (valid) {
                        ArrayList<Couple<Integer, Integer>> outterHexagons = new ArrayList<>();
                        Integer[] occurence = new Integer[pattern.getNbNodes()];

                        for (int i = 0; i < coords.length; i++) {
                            Couple<Integer, Integer> coord = coords[i];
                            int indexOutterHexagon = diameter * diameter;
                            if (coord.getX() >= 0 && coord.getX() < diameter && coord.getY() >= 0
                                    && coord.getY() < diameter) {
                                occurence[i] = hexagonIndicesMatrix[coord.getY()][coord.getX()];

                                if (hexagonIndicesMatrix[coord.getY()][coord.getX()] == -1 && !outterHexagons.contains(coord)) {
                                    outterHexagons.add(coord);
                                    outterHexagonsIndexes.add(indexOutterHexagon);
                                    indexOutterHexagon++;
                                }
                            } else {
                                occurence[i] = -1;

                                if (!outterHexagons.contains(coord)) {
                                    outterHexagons.add(coord);
                                    outterHexagonsIndexes.add(indexOutterHexagon);
                                    indexOutterHexagon++;
                                }
                            }
                        }

                        ArrayList<Integer> present = new ArrayList<>();
                        ArrayList<Integer> absent = new ArrayList<>();
                        ArrayList<Integer> unknown = new ArrayList<>();
                        ArrayList<Integer> outter = new ArrayList<>();

                        for (int i = 0; i < pattern.getNbNodes(); i++) {
                            if (pattern.getLabel(i) == PatternLabel.NEUTRAL) {
                                Couple<Integer, Integer> coord = coords[i];

                                if (coord.getX() >= 0 && coord.getX() < diameter && coord.getY() >= 0
                                        && coord.getY() < diameter) {
                                    if (hexagonIndicesMatrix[coord.getY()][coord.getX()] == -1) {

                                        int index = findIndex(outterHexagons, coord);
                                        outter.add(outterHexagonsIndexes.get(index));
                                    } else {
                                        unknown.add(hexagonIndicesMatrix[coord.getY()][coord.getX()]);
                                    }
                                } else {
                                    int index = findIndex(outterHexagons, coord);
                                    outter.add(outterHexagonsIndexes.get(index));
                                }
                            } else if (pattern.getLabel(i) == PatternLabel.POSITIVE) {
                                Couple<Integer, Integer> coord = coords[i];
                                present.add(hexagonIndicesMatrix[coord.getY()][coord.getX()]);
                            } else if (pattern.getLabel(i) == PatternLabel.NEGATIVE) {
                                Couple<Integer, Integer> coord = coords[i];

                                if (coord.getX() >= 0 && coord.getX() < diameter && coord.getY() >= 0
                                        && coord.getY() < diameter) {
                                    if (hexagonIndicesMatrix[coord.getY()][coord.getX()] == -1) {

                                        int index = findIndex(outterHexagons, coord);
                                        outter.add(outterHexagonsIndexes.get(index));
                                    } else {
                                        absent.add(hexagonIndicesMatrix[coord.getY()][coord.getX()]);
                                    }
                                } else {
                                    int index = findIndex(outterHexagons, coord);
                                    outter.add(outterHexagonsIndexes.get(index));
                                }
                            }
                        }

                        fragmentOccurences.addOccurence(occurence);
                        fragmentOccurences.addCoordinate(coords);
                        fragmentOccurences.addOutterHexagons(outter);
                        fragmentOccurences.addPresentHexagons(present);
                        fragmentOccurences.addAbsentHexagons(absent);
                        fragmentOccurences.addUnknownHexagons(unknown);
                    }
                }
            }
        }

        return fragmentOccurences;
    }

    public ArrayList<ArrayList<Integer>> getNoGoods() {
        return nogoods;
    }

    public BoolVar getNbHexagonsReified(int index) {
        return nbHexagonsReifies[index];
    }

    public void setNbHexagonsReified(int index, BoolVar value) {
        nbHexagonsReifies[index] = value;
    }

    public Model getChocoModel() {
        return chocoModel;
    }

    public SimpleIntegerProperty getNbTotalSolutions() {
        return nbTotalSolutions;
    }

    public ModelPropertySet getModelPropertySet() {
        return this.modelPropertySet;
    }

    public static SolverPropertySet getSolverPropertySet() {
        return GeneralModel.solverPropertySet;
    }

    public static void buildSolverPropertySet(ArrayList<HBoxCriterion> hBoxesSolverCriterions) {
        solverPropertySet.clearPropertyExpressions();
        for (HBoxCriterion box : hBoxesSolverCriterions) {
            ((HBoxSolverCriterion) box).addPropertyExpression(solverPropertySet);
        }
    }

    public void increaseDegree(String variableName) {
        variablesDegrees.putIfAbsent(variableName, 0);
        int curentDegree = variablesDegrees.get(variableName);
        variablesDegrees.remove(variableName);
        variablesDegrees.put(variableName, curentDegree + 1);
    }

    public void displayDegrees() {
        for (Map.Entry<String, Integer> entry : variablesDegrees.entrySet()) {
            String key = entry.getKey();
            Object value = entry.getValue();

            System.out.println("d(" + key + ") = " + value);
        }
    }

    /*
     * Getters & Setters
     */

    public int getNbCrowns() {
        return nbCrowns;
    }

    public int getDiameter() {
        return diameter;
    }

    public Model getProblem() {
        return chocoModel;
    }

    public int[][] getHexagonIndicesMatrix() {
        return hexagonIndicesMatrix;
    }

    public int getHexagonIndex(int lineIndex, int columnIndex){
        return hexagonIndicesMatrix[lineIndex][columnIndex];
    }

    public int[][] getSideSharing() {
        return sideSharing;
    }

    public boolean sharesSide(int i, int j) {
        return sideSharing[i][j] == 1;
    }

    public BoolVar[] getBenzenoidVerticesBVArray() {
        return benzenoidVerticesBVArray;
    }

    public BoolVar getBenzenoidVerticesBVArray(int index) {
        return benzenoidVerticesBVArray[index];
    }

    public UndirectedGraphVar getGraphVar() {
        return benzenoidGraphVar;
    }

    public BoolVar[][] getBenzenoidEdges() {
        return benzenoidEdges;
    }

    public void setInTestMode(boolean inTestMode) {
        isInTestMode = inTestMode;
    }
    // ADDED
    // Dans GeneralModel.java

    public void postTransformationConstraints() {
        if (transformationActiveVars == null || transformationActiveVars.length == 0) {
            // Si aucune transformation n'est possible ou configurée,
            // et que l'utilisateur a demandé des C5/C7, cela devrait être une contradiction.
            if (nbPentagonsVar != null && nbPentagonsVar.getUB() > 0) {
                getProblem().arithm(nbPentagonsVar, "=", 0).post(); // Force 0 C5
            }
            if (nbHeptagonsVar != null && nbHeptagonsVar.getUB() > 0) {
                getProblem().arithm(nbHeptagonsVar, "=", 0).post(); // Force 0 C7
            }
            return;
        }

        postTransformationAtomExistsConstraints();
        postTransformationNonOverlapConstraints();
        // Les contraintes sur nbPentagonsVar et nbHeptagonsVar sont postées
        // par PentagonNumberConstraint et HeptagonNumberConstraint
    }



    private void postTransformationNonOverlapConstraints() {
        for (int k = 0; k < potentialTransformationSites.size(); k++) {
            for (int j = k + 1; j < potentialTransformationSites.size(); j++) {
                Couple<Integer, Integer> site_k = potentialTransformationSites.get(k);
                Couple<Integer, Integer> site_j = potentialTransformationSites.get(j);

                if (sitesAreConflicting(site_k, site_j)) {
                    // Si les sites k et j sont en conflit, ils ne peuvent pas être actifs en même temps.
                    // transformationActiveVars[k] + transformationActiveVars[j] <= 1
                    getProblem().arithm(transformationActiveVars[k], "+", transformationActiveVars[j], "<=", 1).post();
                }
            }
        }
    }

    private boolean sitesAreConflicting(Couple<Integer, Integer> site_k, Couple<Integer, Integer> site_j) {
        // site_k = (h_a, h_b), site_j = (h_c, h_d)
        // Conflit simple: si les sites partagent un hexagone de base.
        if (site_k.getX().equals(site_j.getX()) || site_k.getX().equals(site_j.getY()) ||
                site_k.getY().equals(site_j.getX()) || site_k.getY().equals(site_j.getY())) {
            return true;
        }

        // TODO: Logique plus avancée si nécessaire :
        // Vérifier si un hexagone du site_k est adjacent à un hexagone du site_j
        // (différent de l'adjacence qui définit le site lui-même).
        // Par exemple, si h_a est adjacent à h_c (et (h_a,h_c) n'est pas site_j, etc.)
        // La définition de "conflit" ici est cruciale et dépend de la géométrie de la transformation.
        // Pour l'instant, on se limite au partage d'un hexagone de base.
        return false;
    }
}
