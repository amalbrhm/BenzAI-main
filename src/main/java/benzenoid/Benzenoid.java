package benzenoid;

import classifier.Irregularity;
import generator.patterns.Pattern;
import generator.patterns.PatternLabel;
import generator.properties.Property;
import generator.properties.model.ModelProperty;
import generator.properties.model.ModelPropertySet;
import benzenoid.sort.MoleculeComparator;
import benzenoid.sort.NbHexagonsComparator;
import parsers.GraphCoordFileBuilder;
import parsers.GraphFileBuilder;
import parsers.GraphParser;
import solution.ClarCoverSolution;
import solveur.Aromaticity;
import utils.Couple;
import utils.HexNeighborhood;
import utils.Interval;
import utils.RelativeMatrix;

import java.io.File;
import java.io.IOException;
import java.util.*;
import benzenoid.CycleType; // Ajout de l'import pour le nouveau enum

public class Benzenoid implements Comparable<Benzenoid> {

	// AJOUT: Champ pour stocker les types de cycle des hexagones
	private CycleType[] hexagonCycleTypes;

	private MoleculeComparator comparator;

	private final int nbCarbons;
	private final int nbBonds;
	private final int nbHexagons;
	private int benzdbId;
	private int nbStraightEdges;
	private int maxIndex;
	private ArrayList<ArrayList<Integer>> edgeLists;
	private final int[][] edgeMatrix;
	private ArrayList<String> edgesString;
	private final ArrayList<String> hexagonsString;
	private final Node[] nodesCoordinates;
	private final RelativeMatrix matrixCoordinates;
	private final int[][] hexagons;
	private int[][] dualGraph;
	private int[] degrees;

	private ArrayList<ArrayList<Integer>> hexagonsVertices;

	private int nbHydrogens;

	private String name;
	private String inchi;

	private String description;

	private ArrayList<Integer> verticesSolutions;

	private Couple<Integer, Integer>[] hexagonsCoords;

	private int[][] fixedBonds;
	private int[] fixedCircles;

	private ArrayList<String> names;

	private int nbCrowns = -1;

	private final BenzenoidDatabaseInformation databaseInformation;
	private boolean databaseCheck;    // true if the database has been check, false otherwise

	private final BenzenoidComputableInformations computableInformations;
	private int clarNumber;
	private double homo;
	private double lumo;
	private double moment;


	/**
	 * Constructors
	 */

	/**
	 * Constructors
	 */

	public Benzenoid(int nbNodes, int nbEdges, int nbHexagons, int[][] inputHexagons, Node[] nodesRefs,
					 int[][] edgeMatrix, RelativeMatrix coordinates) {

		comparator = new NbHexagonsComparator();

		this.nbCarbons = nbNodes;
		this.nbBonds = nbEdges;
		this.nbHexagons = nbHexagons;
		// MODIFIÉ: Assigner inputHexagons à this.hexagons. La taille des faces sera utilisée plus tard.
		this.hexagons = inputHexagons;
		this.nodesCoordinates = nodesRefs;
		this.edgeMatrix = edgeMatrix;
		this.matrixCoordinates = coordinates;
		this.inchi = "";
		this.benzdbId = -1;
		this.clarNumber = -1;
		this.homo = 1000000;
		this.lumo = 1000000;
		this.moment = 1000000;

		// AJOUT/MODIFIÉ: Initialisation de hexagonCycleTypes
		if (nbHexagons > 0) {
			this.hexagonCycleTypes = new CycleType[nbHexagons];
			for (int i = 0; i < nbHexagons; i++) {
				this.hexagonCycleTypes[i] = CycleType.UNKNOWN; // Sera défini par determineAndSetCycleTypes
			}
		} else {
			this.hexagonCycleTypes = new CycleType[0];
		}

		hexagonsString = new ArrayList<>();
		// MODIFIÉ: Construction de hexagonsString basée sur this.hexagons (qui est inputHexagons ici)
		if (this.hexagons != null && nodesRefs != null) {
			for (int[] hexagon_nodes : this.hexagons) {
				if (hexagon_nodes != null) {
					StringBuilder builder = new StringBuilder();
					builder.append("h ");
					for (int u : hexagon_nodes) {
						if (u >= 0 && u < nodesRefs.length && nodesRefs[u] != null) {
							Node node = nodesRefs[u];
							builder.append(node.getX()).append("_").append(node.getY()).append(" ");
						}
					}
					this.hexagonsString.add(builder.toString());
				}
			}
		}

		initHexagons(); // initHexagons va utiliser this.hexagonsString pour reconstruire this.hexagons
		// avec les tailles correctes pour this.hexagons[i] si ce n'était pas déjà le cas
		// ou pour peupler hexagonsVertices
		determineAndSetCycleTypes(); // AJOUTÉ: Appel pour définir les types de cycles

		computeDualGraph();
		computeDegrees();
		buildHexagonsCoords2();

		edgeLists = new ArrayList<ArrayList<Integer>> ();

		databaseInformation = new BenzenoidDatabaseInformation(this);
		databaseCheck = false;
		computableInformations = new BenzenoidComputableInformations(this);
	}

	public Benzenoid(int nbNodes, int nbEdges, int nbHexagons, ArrayList<ArrayList<Integer>> edgeLists,
					 int[][] edgeMatrix, ArrayList<String> edgesString, ArrayList<String> hexagonsStringInput,
					 Node[] nodesRefs, RelativeMatrix coords, int maxIndex) {

		this.nbCarbons = nbNodes;
		this.nbBonds = nbEdges;
		this.nbHexagons = nbHexagons;
		this.edgeLists = edgeLists;
		this.edgeMatrix = edgeMatrix;
		this.edgesString = edgesString;
		this.hexagonsString = hexagonsStringInput;
		this.nodesCoordinates = nodesRefs;
		this.matrixCoordinates = coords;
		this.inchi = "";
		this.benzdbId = -1;

		this.maxIndex = maxIndex;

		// AJOUT/MODIFIÉ: Initialisation de this.hexagons et this.hexagonCycleTypes
		if (nbHexagons > 0) {
			this.hexagonCycleTypes = new CycleType[nbHexagons];
			for (int i = 0; i < nbHexagons; i++) {
				this.hexagonCycleTypes[i] = CycleType.UNKNOWN;
			}
			// MODIFIÉ: Initialiser la première dimension de this.hexagons.
			// La deuxième dimension (taille de chaque face) sera définie dans initHexagons.
			this.hexagons = new int[nbHexagons][];
		} else {
			this.hexagonCycleTypes = new CycleType[0];
			this.hexagons = new int[0][];
		}

		initHexagons(); // Peuple this.hexagons[][] (y compris la 2ème dimension) à partir de this.hexagonsString
		determineAndSetCycleTypes(); // AJOUTÉ: Appel après initHexagons

		nbStraightEdges = 0;
		if (this.edgeMatrix != null && nodesRefs != null) {
			for (int i = 0; i < this.edgeMatrix.length; i++) {
				if (i >= nodesRefs.length || nodesRefs[i] == null) continue;
				Node u1 = nodesRefs[i];

				for (int j = (i + 1); j < this.edgeMatrix[i].length; j++) {
					if (j >= nodesRefs.length || nodesRefs[j] == null) continue;
					Node u2 = nodesRefs[j];

					if (this.edgeMatrix[i][j] == 1) {
						if (u1.getX() == u2.getX())
							nbStraightEdges++;
					}
				}
			}
		}

		computeDualGraph();
		computeDegrees();
		buildHexagonsCoords2();

		databaseInformation = new BenzenoidDatabaseInformation(this);
		databaseCheck = false;
		computableInformations = new BenzenoidComputableInformations(this);
	}
	// AJOUTÉ: Nouvelle méthode pour déterminer et définir les types de cycle
	private void determineAndSetCycleTypes() {
		if (this.hexagons == null || this.hexagonCycleTypes == null || this.nbHexagons == 0) {
			return;
		}
		for (int i = 0; i < this.nbHexagons; i++) {
			if (this.hexagons[i] != null) {

				int numberOfNodesInFace = this.hexagons[i].length;

				switch (numberOfNodesInFace) {
					case 5:
						this.setHexagonCycleType(i, CycleType.C5);
						break;
					case 6:
						this.setHexagonCycleType(i, CycleType.C6);
						break;
					case 7:
						this.setHexagonCycleType(i, CycleType.C7);
						break;
					default:
						this.setHexagonCycleType(i, CycleType.UNKNOWN);
						break;
				}
			} else {
				this.setHexagonCycleType(i, CycleType.UNKNOWN);
			}
		}
	}
	/**
	 * Getters and setters
	 */

	public int[][] getDualGraph() {
		return dualGraph;
	}

	public int getNbCarbons() {
		return nbCarbons;
	}

	public int getNbBonds() {
		return nbBonds;
	}

	public int getNbHexagons() {
		return nbHexagons;
	}

	public int getMaxIndex() {
		return maxIndex;
	}


	public ArrayList<ArrayList<Integer>> getEdgeLists() {
		return edgeLists;
	}

	public int[][] getEdgeMatrix() {
		return edgeMatrix;
	}

	public ArrayList<String> getEdgesString() {
		return edgesString;
	}

	public ArrayList<String> getHexagonsString() {
		return hexagonsString;
	}

	public Node getNodeRef(int index) {
		return nodesCoordinates[index];
	}

	public RelativeMatrix getMatrixCoordinates() {
		return matrixCoordinates;
	}

	public Node[] getNodesCoordinates() {
		return nodesCoordinates;
	}

	public ArrayList<ArrayList<Integer>> getHexagonsVertices() {
		return hexagonsVertices;
	}

	public int getNbStraightEdges() {
		return nbStraightEdges;
	}

	public int[][] getHexagons() {
		return hexagons;
	}

	public int degree(int u) {
		return degrees[u];
	}

	public boolean hasCheckedDatabase() {
		return databaseCheck;
	}

	public void performCheckDatabase() {
		databaseCheck = true;
	}

	/**
	 * Class's methods
	 */

	private void computeDegrees() {

		degrees = new int[nbCarbons];

		for (int i = 0; i < nbCarbons; i++) {

			int degree = 0;
			for (int j = 0; j < nbCarbons; j++) {

				if (edgeMatrix[i][j] == 1)
					degree++;
			}

			degrees[i] = degree;
		}
	}

	private void computeDualGraph() {

		dualGraph = new int[nbHexagons][6];

		for (int i = 0; i < nbHexagons; i++)
			Arrays.fill(dualGraph[i], -1);

		ArrayList<Integer> candidats = new ArrayList<>();
		candidats.add(0);

		int index = 0;

		while (index < nbHexagons) {

			int candidat = candidats.get(index);
			int[] candidatHexagon = hexagons[candidat];

			for (int i = 0; i < candidatHexagon.length; i++) {

				int u = candidatHexagon[i];
				int v = candidatHexagon[(i + 1) % 6];

				// System.out.print(""); // Original line, can be kept or removed

				for (int j = 0; j < nbHexagons; j++) {
					if (j != candidat) {

						int contains = 0;
						for (int k = 0; k < 6; k++) {
							if (hexagons[j][k] == u || hexagons[j][k] == v)
								contains++;
						}

						if (contains == 2) {

							dualGraph[candidat][i] = j;

							if (!candidats.contains(j))
								candidats.add(j);

							break;
						}
					}
				}

			}
			index++;
		}
	}


	private int findHexagon(int u, int v) {

		for (int i = 0; i < nbHexagons; i++) {
			int[] hexagon = hexagons[i];

			if (hexagon[4] == u && hexagon[5] == v)
				return i;
		}

		return -1;

	}

	private ArrayList<Integer> findHexagons(int hexagon, Interval interval) {

		ArrayList<Integer> hexagons = new ArrayList<>();
		int size = interval.size() / 2;

		hexagons.add(hexagon);

		int newHexagon = hexagon;

		for (int i = 0; i < size; i++) {

			newHexagon = dualGraph[newHexagon][1];
			hexagons.add(newHexagon);
		}

		return hexagons;

	}

	public ArrayList<Integer> getAllHexagonsOfIntervals(ArrayList<Interval> intervals) {

		ArrayList<Integer> hexagons = new ArrayList<>();

		for (Interval interval : intervals) {

			int hexagon = findHexagon(interval.x1(), interval.y1());
			hexagons.addAll(findHexagons(hexagon, interval));
		}

		return hexagons;
	}
	// Modifié
	private void initHexagons() {

		hexagonsVertices = new ArrayList<>();

		for (int i = 0; i < nbCarbons; i++) // MODIFIÉ : S'assurer que nbCarbons est utilisé pour la taille ou ajouter une vérification
			hexagonsVertices.add(new ArrayList<>());

		if (hexagonsString == null || matrixCoordinates == null || this.hexagons == null) { // AJOUT: Protection
			return;
		}

		for (int i = 0; i < nbHexagons; i++) {
			if (i >= hexagonsString.size()) { // AJOUT: Protection
				if (i < this.hexagons.length) this.hexagons[i] = new int[0];
				continue;
			}
			String hexagonStrValue = hexagonsString.get(i); // Renommé 'hexagon' en 'hexagonStrValue' pour éviter la confusion
			String[] sHexagon = hexagonStrValue.split(" ");

			// AJOUTÉ/MODIFIÉ : Allouer this.hexagons[i] ici avec la bonne taille
			if (i < this.hexagons.length) { // S'assurer que l'index i est valide
				this.hexagons[i] = new int[sHexagon.length - 1];
			} else {
				continue;
			}

			for (int j = 1; j < sHexagon.length; j++) {
				String[] sVertex = sHexagon[j].split(java.util.regex.Pattern.quote("_"));
				if (sVertex.length < 2) { // AJOUT: Protection
					if (j-1 < this.hexagons[i].length) this.hexagons[i][j-1] = -1;
					continue;
				}
				int x = Integer.parseInt(sVertex[0]);
				int y = Integer.parseInt(sVertex[1]);
				int nodeIndexInMolecule = matrixCoordinates.get(x, y);

				if (nodeIndexInMolecule != -1) {
					if (j-1 < this.hexagons[i].length) {
						this.hexagons[i][j - 1] = nodeIndexInMolecule;
					}
					// AJOUT: Protection pour hexagonsVertices
					if (nodeIndexInMolecule >= 0 && nodeIndexInMolecule < hexagonsVertices.size()) {
						hexagonsVertices.get(nodeIndexInMolecule).add(i);
					}
				} else {
					if (j-1 < this.hexagons[i].length) this.hexagons[i][j-1] = -1;
				}
			}
		}
	}

	public int getNbHydrogens() {

		if (nbHydrogens == 0) {

			for (int i = 0; i < nbCarbons; i++) {

				int degree = 0;
				for (int j = 0; j < nbCarbons; j++) {

					if (edgeMatrix[i][j] == 1)
						degree++;
				}

				if (degree == 2)
					nbHydrogens++;
			}
		}

		return nbHydrogens;
	}

	@SuppressWarnings("unchecked")
	@Override
	public String toString() {

		if (name == null) {

			int nbCrowns = (int) Math.floor((((double) nbHexagons + 1) / 2.0) + 1.0);

			if (nbHexagons % 2 == 1)
				nbCrowns--;

			int diameter = (2 * nbCrowns) - 1;

			int[][] coordsMatrix = new int[diameter][diameter];

			/*
			 * Building coords matrix
			 */

			for (int i = 0; i < diameter; i++) {
				Arrays.fill(coordsMatrix[i], -1);
			}

			//for (int i = 0; i < diameter; i++) // Ligne redondante
			//	Arrays.fill(coordsMatrix[i], -1);

			int index = 0;
			int m = (diameter - 1) / 2;

			int shift = diameter - nbCrowns;

			for (int i = 0; i < m; i++) {

				for (int j = 0; j < diameter - shift; j++) {
					coordsMatrix[i][j] = index;
					index++;
				}

				for (int j = diameter - shift; j < diameter; j++)
					index++;

				shift--;
			}

			for (int j = 0; j < diameter; j++) {
				coordsMatrix[m][j] = index;
				index++;
			}

			shift = 1;

			for (int i = m + 1; i < diameter; i++) {

				for (int j = 0; j < shift; j++)
					index++;

				for (int j = shift; j < diameter; j++) {
					coordsMatrix[i][j] = index;
					index++;
				}

				shift++;
			}

			Couple<Integer, Integer>[] localHexagonsCoords = new Couple[nbHexagons]; // Renommé pour éviter la confusion avec le champ de classe

			ArrayList<Integer> candidats = new ArrayList<>();
			candidats.add(0);

			localHexagonsCoords[0] = new Couple<>(0, 0); // Utilisation de la variable locale

			int[] checkedHexagons = new int[nbHexagons];
			checkedHexagons[0] = 1;

			while (!candidats.isEmpty()) {
				int candidat = candidats.get(0);
				for (HexNeighborhood neighbor : HexNeighborhood.values()) {
					int neighborIndex = dualGraph[candidat][neighbor.getIndex()];
					if (neighborIndex != -1) {
						if (checkedHexagons[neighborIndex] == 0) {
							localHexagonsCoords[neighborIndex] = new Couple<>(localHexagonsCoords[candidat].getX() + neighbor.dx(), localHexagonsCoords[candidat].getY() + neighbor.dy());
							candidats.add(neighborIndex);
							checkedHexagons[neighborIndex] = 1;
						}
					}
				}

				candidats.remove(candidats.get(0));
			}

			/*
			 * Trouver une manière de le faire rentrer dans le coronénoïde
			 */

			StringBuilder code = new StringBuilder();

			for (int i = 0; i < diameter; i++) {
				for (int j = 0; j < diameter; j++) {

					int h = coordsMatrix[i][j];

					if (h != -1) {

						boolean ok = true;

						Couple<Integer, Integer>[] newHexagons = new Couple[localHexagonsCoords.length];

						for (int k = 0; k < localHexagonsCoords.length; k++) {

							Couple<Integer, Integer> hexagon = localHexagonsCoords[k]; // Utilisation de la variable locale

							int xh = hexagon.getY() + i;
							int yh = hexagon.getX() + j;

							if (xh >= 0 && xh < diameter && yh >= 0 && yh < diameter && coordsMatrix[xh][yh] != -1) {
								newHexagons[k] = new Couple<>(xh, yh);
							}

							else {
								ok = false;
								break;
							}
						}

						if (ok) {

							for (int k = 0; k < newHexagons.length; k++) {
								code.append(coordsMatrix[newHexagons[k].getX()][newHexagons[k].getY()]);

								if (k < newHexagons.length - 1)
									code.append("-");
							}

							name = code.toString();
							return code.toString();

						}

					}
				}
			}

			name = null;
			return null;
		}

		else
			return name;
	}

	public void setVerticesSolutions(ArrayList<Integer> verticesSolutions) {
		this.verticesSolutions = verticesSolutions;
	}

	public ArrayList<Integer> getVerticesSolutions() {
		return verticesSolutions;
	}

	// AJOUTÉ/MODIFIÉ: Getters et Setters pour hexagonCycleTypes avec protections
	public CycleType getHexagonCycleType(int hexagonIndex) {
		if (hexagonCycleTypes != null && hexagonIndex >= 0 && hexagonIndex < this.hexagonCycleTypes.length ) {
			return this.hexagonCycleTypes[hexagonIndex];
		}
		return CycleType.UNKNOWN;
	}

	public void setHexagonCycleType(int hexagonIndex, CycleType cycleType) {
		if (hexagonCycleTypes != null && hexagonIndex >= 0 && hexagonIndex < this.hexagonCycleTypes.length ) {
			this.hexagonCycleTypes[hexagonIndex] = cycleType;
		}
	}

	public void setDescription(String description) {
		this.description = description;
	}

	public String getDescription() {
		return description;
	}

	public Couple<Integer, Integer>[] getHexagonsCoords() {
		return hexagonsCoords;
	}

	public int[] getDegrees() {
		return degrees;
	}

	public Irregularity getIrregularity() {
		return computableInformations.computeParameterOfIrregularity();
	}

	public Optional<Aromaticity> getAromaticity() {
		return computableInformations.getAromaticity();
	}

	public double getNbKekuleStructures() {
		return computableInformations.getNbKekuleStructures();
	}

	public boolean isAromaticityComputed() {
		return computableInformations.isAromaticityComputed();
	}

	public boolean isClarCoverComputed() {
		return computableInformations.isClarCoverComputed();
	}

	public boolean edgeExists(int i, int j) {
		return edgeMatrix[i][j] == 1;
	}

	public ArrayList<Couple<Integer, Integer>> getBoundsInvolved(int carbon) {

		ArrayList<Couple<Integer, Integer>> bounds = new ArrayList<>();

		for (int i = 0; i < nbCarbons; i++) {
			if (edgeExists(carbon, i))
				bounds.add(new Couple<>(carbon, i));
		}

		return bounds;

	}

	public ArrayList<Integer> getHexagonsInvolved(int carbon) {

		ArrayList<Integer> hexagonsList = new ArrayList<>(); // Renommé pour éviter la confusion

		for (int i = 0; i < nbHexagons; i++) {
			int[] hexagon = getHexagons()[i];
			for (Integer u : hexagon) {
				if (u == carbon) {
					hexagonsList.add(i);
					break;
				}
			}
		}

		return hexagonsList;
	}

	public ArrayList<Integer> getHexagonsInvolved(int carbon1, int carbon2) {

		ArrayList<Integer> hexagonsList = new ArrayList<>(); // Renommé pour éviter la confusion

		for (int i = 0; i < nbHexagons; i++) {

			int[] hexagon = getHexagons()[i];
			boolean containsCarbon1 = false;
			boolean containsCarbon2 = false;

			for (Integer u : hexagon) {
				if (u == carbon1)
					containsCarbon1 = true;
				if (u == carbon2)
					containsCarbon2 = true;

				if (containsCarbon1 && containsCarbon2) { // Correction: ajouter l'hexagon APRES avoir vérifié les deux
					hexagonsList.add(i);
					break; // Sortir de la boucle interne une fois l'hexagon trouvé
				}
			}
		}
		return hexagonsList;
	}

	public int[] getHexagon(int hexagon) {
		return hexagons[hexagon];
	}

	public RBO getRBO() {
		return computableInformations.getRingBondOrder();
	}

	public void setComparator(MoleculeComparator comparator) {
		this.comparator = comparator;
	}

	@Override
	public int compareTo(Benzenoid arg0) {
		return comparator.compare(this, arg0);
	}

	public void setClarCoverSolution(ClarCoverSolution clarCoverSolution) {
		computableInformations.setClarCoverSolution(clarCoverSolution);
	}

	public ClarCoverSolution getClarCoverSolution() {
		return computableInformations.getClarCoverSolution();
	}

	public int[] resonanceEnergyClar() {
		return computableInformations.clarResonanceEnergy();
	}

	private int[][] buildCoordsMatrix(int nbCrowns, int diameter) {

		int[][] coordsMatrix = new int[diameter][diameter];
		for (int i = 0; i < diameter; i++) {
			Arrays.fill(coordsMatrix[i], -1);
		}

		int index = 0;
		int m = (diameter - 1) / 2;

		int shift = diameter - nbCrowns;

		for (int i = 0; i < m; i++) {

			for (int j = 0; j < diameter - shift; j++) {
				coordsMatrix[i][j] = index;
				index++;
			}

			for (int j = diameter - shift; j < diameter; j++)
				index++;

			shift--;
		}

		for (int j = 0; j < diameter; j++) {
			coordsMatrix[m][j] = index;
			index++;
		}

		shift = 1;

		for (int i = m + 1; i < diameter; i++) {

			for (int j = 0; j < shift; j++)
				index++;

			for (int j = shift; j < diameter; j++) {
				coordsMatrix[i][j] = index;
				index++;
			}

			shift++;
		}

		return coordsMatrix;
	}

	private boolean touchBorder(int xShift, int yShift, ArrayList<Integer> topBorder, ArrayList<Integer> leftBorder,
								int[][] coordsMatrix) {

		boolean topContains = false;
		boolean leftContains = false;

		for (int i = 0; i < nbHexagons; i++) {
			Couple<Integer, Integer> coord = hexagonsCoords[i];

			int x = coord.getX() + xShift;
			int y = coord.getY() + yShift;

			if (x >= 0 && x < coordsMatrix.length && y >= 0 && y < coordsMatrix.length) {

				if (topBorder.contains(coordsMatrix[x][y]))
					topContains = true;

				if (leftBorder.contains(coordsMatrix[x][y]))
					leftContains = true;

				if (topContains && leftContains)
					return true;
			}
		}

		return false;
	}

	private boolean touchBorder(ArrayList<Integer> hexagonsInput, ArrayList<Integer> topBorder, // Renommé pour éviter la confusion
								ArrayList<Integer> leftBorder) {

		boolean touchTop = false;
		boolean touchLeft = false;

		for (Integer hexagon : hexagonsInput) { // Utilisation du paramètre renommé
			if (topBorder.contains(hexagon))
				touchTop = true;
			if (leftBorder.contains(hexagon))
				touchLeft = true;
		}

		return touchTop && touchLeft;
	}

	private boolean areNeighbors(int hexagon1, int hexagon2) {

		for (int i = 0; i < 6; i++) {
			if (dualGraph[hexagon1][i] == hexagon2)
				return true;
		}

		return false;
	}

	public Pattern convertToPattern(int xShift, int yShift) {

		@SuppressWarnings("unchecked")
		Couple<Integer, Integer>[] shiftCoords = new Couple[nbHexagons];

		for (int i = 0; i < nbHexagons; i++)
			shiftCoords[i] = new Couple<>(hexagonsCoords[i].getX() + xShift, hexagonsCoords[i].getY() + yShift);

		/*
		 * Nodes
		 */

		Node[] patternNodes = new Node[nbHexagons]; // Renommé pour éviter confusion avec le champ de classe

		for (int i = 0; i < nbHexagons; i++) {
			Couple<Integer, Integer> couple = shiftCoords[i];
			patternNodes[i] = new Node(couple.getX(), couple.getY(), i); // Utilisation de la variable locale
		}

		/*
		 * Matrix
		 */

		int[][] matrix = new int[nbHexagons][nbHexagons];
		int[][] neighbors = new int[nbHexagons][6];

		for (int i = 0; i < nbHexagons; i++)
			Arrays.fill(neighbors[i], -1);

		for (int i = 0; i < nbHexagons; i++) {

			Node n1 = patternNodes[i]; // Utilisation de la variable locale

			for (int j = (i + 1); j < nbHexagons; j++) {

				Node n2 = patternNodes[j]; // Utilisation de la variable locale

				if (areNeighbors(i, j)) {

					// Setting matrix
					matrix[i][j] = 1;
					matrix[j][i] = 1;

					// Setting neighbors
					int x1 = n1.getX();
					int y1 = n1.getY();
					int x2 = n2.getX();
					int y2 = n2.getY();

					if (x2 == x1 && y2 == y1 - 1) {
						neighbors[i][0] = j;
						neighbors[j][3] = i;
					}

					else if (x2 == x1 + 1 && y2 == y1) {
						neighbors[i][1] = j;
						neighbors[j][4] = i;
					}

					else if (x2 == x1 + 1 && y2 == y1 + 1) {
						neighbors[i][2] = j;
						neighbors[j][5] = i;
					}

					else if (x2 == x1 && y2 == y1 + 1) {
						neighbors[i][3] = j;
						neighbors[j][0] = i;
					}

					else if (x2 == x1 - 1 && y2 == y1) {
						neighbors[i][4] = j;
						neighbors[j][1] = i;
					}

					else if (x2 == x1 - 1 && y2 == y1 - 1) {
						neighbors[i][5] = j;
						neighbors[j][2] = i;
					}
				}

			}
		}

		/*
		 * Label
		 */

		PatternLabel[] labels = new PatternLabel[nbHexagons];

		for (int i = 0; i < nbHexagons; i++)
			labels[i] = PatternLabel.POSITIVE;

		return new Pattern(matrix, labels, patternNodes, null, neighbors, 0); // Utilisation de la variable locale
	}

	private ArrayList<String> translations(Pattern pattern, int diameter, int[][] coordsMatrix,
										   ArrayList<Integer> topBorder, ArrayList<Integer> leftBorder) {

		ArrayList<String> localNames = new ArrayList<>(); // Renommé pour éviter confusion

		int xShiftMax = Math.max(Math.abs(pattern.xMax() - pattern.xMin()), diameter);
		int yShiftMax = Math.max(Math.abs(pattern.yMax() - pattern.yMin()), diameter);

		for (int xShift = 0; xShift < xShiftMax; xShift++) {
			for (int yShift = 0; yShift < yShiftMax; yShift++) {

				Node[] initialCoords = pattern.getNodesRefs();
				Node[] shiftedCoords = new Node[initialCoords.length];

				boolean ok = true;

				for (int i = 0; i < shiftedCoords.length; i++) {
					Node node = initialCoords[i];
					Node newNode = new Node(node.getX() + xShift, node.getY() + yShift, i);
					shiftedCoords[i] = newNode;

					int x = newNode.getX();
					int y = newNode.getY();

					if (x < 0 || x >= diameter || y < 0 || y >= diameter || coordsMatrix[x][y] == -1) {
						ok = false;
						break;
					}
				}

				if (ok) {

					ArrayList<Integer> localHexagons = new ArrayList<>(); // Renommé
					for (Node node : shiftedCoords) {
						int hexagon = coordsMatrix[node.getX()][node.getY()];
						localHexagons.add(hexagon); // Utilisation de la variable locale
					}

					if (touchBorder(localHexagons, topBorder, leftBorder)) { // Utilisation de la variable locale
						Collections.sort(localHexagons); // Utilisation de la variable locale
						StringBuilder nameBuilder = new StringBuilder(); // Renommé
						for (int i = 0; i < localHexagons.size(); i++) { // Utilisation de la variable locale
							nameBuilder.append(localHexagons.get(i)); // Utilisation de la variable locale
							if (i < localHexagons.size() - 1) // Utilisation de la variable locale
								nameBuilder.append("-");
						}
						localNames.add(nameBuilder.toString()); // Utilisation de la variable locale
					}

				}
			}
		}

		return localNames; // Utilisation de la variable locale
	}

	public ArrayList<String> getNames() {

		if (names != null)
			return names;

		/*
		 * Creating coronenoid matrix
		 */

		names = new ArrayList<>();

		int localNbCrowns = (int) Math.floor((((double) nbHexagons + 1) / 2.0) + 1.0); // Renommé pour éviter la confusion

		if (nbHexagons % 2 == 1)
			localNbCrowns--; // Renommé

		int diameter = (2 * localNbCrowns) - 1; // Renommé

		int[][] coordsMatrix = buildCoordsMatrix(localNbCrowns, diameter); // Renommé

		ArrayList<Integer> topBorder = new ArrayList<>();
		ArrayList<Integer> leftBorder = new ArrayList<>();

		for (int i = 0; i < diameter; i++) {
			if (coordsMatrix[0][i] != -1)
				topBorder.add(coordsMatrix[0][i]);
		}

		for (int i = 0; i < diameter; i++) {

			int j = 0;
			while (coordsMatrix[i][j] == -1)
				j++;

			leftBorder.add(coordsMatrix[i][j]);
		}

		int xShift = -1;
		int yShift = -1;

		boolean touchBorderFound = false; // Renommé
		for (int i = 0; i < diameter; i++) {
			for (int j = 0; j < diameter; j++) {
				if (touchBorder(i, j, topBorder, leftBorder, coordsMatrix)) {
					touchBorderFound = true; // Renommé
					xShift = i;
					yShift = j;
					break;
				}
			}
			if (touchBorderFound) // Renommé
				break;
		}

		Pattern pattern = convertToPattern(xShift, yShift);
		ArrayList<Pattern> rotations = pattern.computeRotations();

		for (Pattern f : rotations) {
			ArrayList<String> rotationsNames = translations(f, diameter, coordsMatrix, topBorder, leftBorder);
			for (String nameStr : rotationsNames) { // Renommé pour éviter confusion
				if (!names.contains(nameStr))
					names.add(nameStr);
			}
		}

		return names;
	}

	@SuppressWarnings("unchecked")
	private void buildHexagonsCoords2() {
		int[][] dualGraph = this.getDualGraph();

		int[] checkedHexagons = new int[this.getNbHexagons()];

		ArrayList<Integer> candidates = new ArrayList<>();
		candidates.add(0);

		hexagonsCoords = new Couple[nbHexagons];

		checkedHexagons[0] = 1;
		hexagonsCoords[0] = new Couple<>(0, 0);

		while (!candidates.isEmpty()) {

			int candidateIndex = candidates.get(0);

			for (HexNeighborhood neighbor : HexNeighborhood.values()) {
				int n = dualGraph[candidateIndex][neighbor.getIndex()];
				if (n != -1) {
					if (checkedHexagons[n] == 0) {
						int x = hexagonsCoords[candidateIndex].getX() + neighbor.dx();
						int y = hexagonsCoords[candidateIndex].getY() + neighbor.dy();
						checkedHexagons[n] = 1;
						hexagonsCoords[n] = new Couple<>(x, y);
						candidates.add(n);
					}
				}
			}

			candidates.remove(candidates.get(0));
		}

	}

	public static ArrayList<Benzenoid> union(ArrayList<Benzenoid> molecules1, ArrayList<Benzenoid> molecules2) {

		ArrayList<Benzenoid> molecules = new ArrayList<>();

		for (Benzenoid molecule : molecules1) {
			if (!molecules.contains(molecule))
				molecules.add(molecule);
		}

		for (Benzenoid molecule : molecules2) {
			if (!molecules.contains(molecule))
				molecules.add(molecule);
		}

		return molecules;
	}

	public static ArrayList<Benzenoid> intersection(ArrayList<Benzenoid> molecules1, ArrayList<Benzenoid> molecules2) {

		ArrayList<Benzenoid> molecules = new ArrayList<>();

		for (Benzenoid molecule : molecules1) {
			if (molecules2.contains(molecule))
				molecules.add(molecule);
		}

		return molecules;
	}

	public static ArrayList<Benzenoid> diff(ArrayList<Benzenoid> molecules1, ArrayList<Benzenoid> molecules2) {

		ArrayList<Benzenoid> molecules = new ArrayList<>();

		for (Benzenoid molecule : molecules1) {
			if (!molecules2.contains(molecule))
				molecules.add(molecule);
		}

		return molecules;
	}

	public void setNbCrowns(int nbCrowns) {
		this.nbCrowns = nbCrowns;
	}

	public int getNbCrowns() {
		return nbCrowns;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		if (obj == null || getClass() != obj.getClass()) return false;
		Benzenoid molecule = (Benzenoid) obj;

		// Ensure names are computed for both objects before comparison
		List<String> thisNames = this.getNames();
		List<String> otherNames = molecule.getNames();

		if (thisNames == null || otherNames == null) return false; // Or throw an exception

		for (String nameStr : otherNames) { // Renommé pour éviter confusion
			if (thisNames.contains(nameStr))
				return true;
		}
		return false;
	}

	public void setClarCoverSolutions(List<ClarCoverSolution> clarCoverSolutions) {
		computableInformations.setClarCoverSolutions(clarCoverSolutions);
	}

	public List<ClarCoverSolution> getClarCoverSolutions() {

		return computableInformations.getClarCoverSolutions();
	}

	public void setFixedBonds(int[][] fixedBonds) {
		this.fixedBonds = fixedBonds;
	}

	public int[][] getFixedBonds() {
		return fixedBonds;
	}

	public int[] getFixedCircles() {
		return fixedCircles;
	}

	public void setFixedCircles(int[] fixedCircles) {
		this.fixedCircles = fixedCircles;
	}

	public int colorShift() {
		int colorShift = 0;
		for (int i = 0; i < nbHexagons; i++) {
			if (dualGraph[i][0] == -1 && dualGraph[i][5] == -1)
				colorShift++;
			if (dualGraph[i][2] == -1 && dualGraph[i][3] == -1)
				colorShift--;
		}
		return colorShift;
	}

	public List<int[][]> getKekuleStructures() {
		return computableInformations.getKekuleStructures();
	}

	public void setKekuleStructures(List<int[][]> kekuleStructures) {
		computableInformations.setKekuleStructures(kekuleStructures);
	}

	/***
	 *
	 */
	public boolean respectPostProcessing(ModelPropertySet modelPropertySet) {
		for(Property property : modelPropertySet)
			if(property.hasExpressions() && !((ModelProperty) property).getChecker().checks(this, (ModelProperty) property))
				return false;
		return true;
	}

	/***
	 *
	 */
	public static Benzenoid buildMolecule(String description, int nbCrowns, int index, ArrayList<Integer> verticesSolution) {

		Benzenoid molecule = null;
		try {
			String graphFilename = "tmp.graph";
			String graphCoordFilename = "tmp.graph_coord";

			buildGraphFile(nbCrowns, verticesSolution, graphFilename);
			convertGraphCoordFileInstance(graphFilename, graphCoordFilename);
			molecule = buildMolecule(description, nbCrowns, index, verticesSolution, graphCoordFilename);
			deleteTmpFiles(graphFilename, graphCoordFilename);
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		return molecule;

	}

	/***
	 *
	 */
	private static Benzenoid buildMolecule(String description, int nbCrowns, int index, ArrayList<Integer> verticesSolution,
										   String graphCoordFilename) {
		Benzenoid molecule = GraphParser.parseUndirectedGraph(graphCoordFilename, null, false);
		molecule.setVerticesSolutions(verticesSolution);
		molecule.setDescription(buildMoleculeDescription(description, index));
		molecule.setNbCrowns(nbCrowns);
		return molecule;
	}
	/***
	 *
	 */
	private static String buildMoleculeDescription(String description, int index) {
		String[] lines = description.split("\n");
		StringBuilder b = new StringBuilder();

		b.append("solution_").append(index).append("\n");
		for (int j = 1; j < lines.length; j++)
			b.append(lines[j]).append("\n");
		return b.toString();
	}
	/***
	 *
	 */
	private static void deleteTmpFiles(String graphFilename, String graphCoordFilename) {
		File file = new File(graphFilename);
		if (file.exists()) file.delete(); // Ajout de la vérification d'existence
		file = new File(graphCoordFilename);
		if (file.exists()) file.delete(); // Ajout de la vérification d'existence
	}

	/***
	 *
	 */
	private static void convertGraphCoordFileInstance(String graphFilename, String graphCoordFilename) {
		GraphCoordFileBuilder graphCoordBuilder = new GraphCoordFileBuilder(graphFilename, graphCoordFilename);
		graphCoordBuilder.convertInstance();
	}

	/***
	 *
	 */
	private static void buildGraphFile(int nbCrowns, ArrayList<Integer> verticesSolution,
									   String graphFilename) throws IOException {
		GraphFileBuilder graphBuilder = new GraphFileBuilder(verticesSolution, graphFilename,
				nbCrowns);
		graphBuilder.buildGraphFile();
	}

	public BenzenoidDatabaseInformation getDatabaseInformation() {
		return databaseInformation;
	}

	public void setInchi (String inchi) {
		this.inchi = inchi;
	}

	public String getInchi () {
		return this.inchi;
	}

	public void setClarNumber (int clarNumber) {
		this.clarNumber = clarNumber;
	}

	public int getClarNumber () {
		if (clarNumber == -1 && computableInformations.getClarCoverSolution() != null) { // Vérifier que la solution Clar existe
			clarNumber = this.computableInformations.getClarCoverSolution().getClarNumber();
		}
		return clarNumber;
	}

	public void setHomo (double homo) {
		this.homo = homo;
	}

	public double getHomo () {
		return this.homo ;
	}

	public void setLumo (double lumo) {
		this.lumo = lumo;
	}

	public double getLumo () {
		return this.lumo ;
	}

	public void setMmoment (double moment) {
		this.moment = moment;
	}

	public double getMoment () {
		return this.moment ;
	}

	public void setBenzdbId(int id) {
		benzdbId = id;
	}

	public int getBenzdbId() {
		return benzdbId;
	}


}