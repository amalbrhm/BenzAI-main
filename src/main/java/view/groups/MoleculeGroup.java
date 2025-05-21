package view.groups;

import javafx.scene.Group;
import javafx.scene.text.Text;
import javafx.scene.text.TextBoundsType;
import benzenoid.Benzenoid;
import benzenoid.CycleType; // Ajout de l'import
import utils.Couple;
import utils.HexNeighborhood;

import java.util.ArrayList;

public class MoleculeGroup extends Group {

	private double width;
	private double height;

	protected Benzenoid molecule;
	protected Couple<Integer, Integer>[] hexagonsCoords;
	protected Couple<Double, Double>[] centersCoords;

	protected Hexagon2[] hexagons;
	protected ArrayList<Text> texts; // Pour les IDs
	protected ArrayList<Text> cycleTypeTexts; // Nouveau : pour les types de cycle

	protected boolean writeText = true;

	protected double xShift;
	protected double yShift;

	public MoleculeGroup(Benzenoid molecule) {
		this.molecule = molecule;
		// Initialisation des nouvelles listes de textes
		texts = new ArrayList<>();
		cycleTypeTexts = new ArrayList<>();
		buildHexagons();
		drawHexagons();
	}

	@SuppressWarnings("unchecked")
	private void buildHexagons() {

		hexagons = new Hexagon2[molecule.getNbHexagons()];
		hexagonsCoords = new Couple[molecule.getNbHexagons()];
		centersCoords = new Couple[molecule.getNbHexagons()];
		int[][] dualGraph = molecule.getDualGraph();

		int[] checkedHexagons = new int[molecule.getNbHexagons()];

		ArrayList<Integer> candidates = new ArrayList<>();
		candidates.add(0);

		checkedHexagons[0] = 1;
		hexagonsCoords[0] = new Couple<>(0, 0);

		centersCoords[0] = new Couple<>(40.0, 40.0);

		while (candidates.size() > 0) {

			int candidateIndex = candidates.get(0);

			for (HexNeighborhood neighbor : HexNeighborhood.values()) {

				int n = dualGraph[candidateIndex][neighbor.getIndex()];
				if (n != -1) {
					if (checkedHexagons[n] == 0) {

						int x = hexagonsCoords[candidateIndex].getX() + neighbor.dx();
						int y = hexagonsCoords[candidateIndex].getY() + neighbor.dy();

						double xCenter = centersCoords[candidateIndex].getX();
						double yCenter = centersCoords[candidateIndex].getY();

						switch (neighbor.getIndex()) {
							case 0:
								xCenter += 26.0;
								yCenter -= 43.5;
								break;
							case 1:
								xCenter += 52.0;
								yCenter += 0.0;
								break;
							case 2:
								xCenter += 26.0;
								yCenter += 43.5;
								break;
							case 3:
								xCenter -= 26.0;
								yCenter += 43.5;
								break;
							case 4:
								xCenter -= 52.0;
								yCenter += 0.0;
								break;
							case 5:
								xCenter -= 26.0;
								yCenter -= 43.5;
						}

						checkedHexagons[n] = 1;
						hexagonsCoords[n] = new Couple<>(x, y);
						centersCoords[n] = new Couple<>(xCenter, yCenter);
						candidates.add(n);
					}
				}
			}

			candidates.remove(candidates.get(0));
		}

		ArrayList<ArrayList<Double>> points = new ArrayList<>();
		for (Couple<Double, Double> centersCoord : centersCoords) {

			points.add(getHexagonPoints(centersCoord.getX(), centersCoord.getY()));
		}

		double xMin = Double.MAX_VALUE;
		double xMax = Double.MIN_VALUE;

		double yMin = Double.MAX_VALUE;
		double yMax = Double.MIN_VALUE;

		for (ArrayList<Double> point : points) {

			for (int j = 0; j < point.size(); j++) {

				double u = point.get(j);

				if (j % 2 == 0) { // x
					if (u < xMin)
						xMin = u;

					if (u > xMax)
						xMax = u;
				} else { // y
					if (u < yMin)
						yMin = u;

					if (u > yMax)
						yMax = u;
				}
			}
		}

		xShift = 0;
		yShift = 0;

		if (xMin < 0)
			xShift = -xMin + 15.0;

		if (yMin < 0)
			yShift = -yMin + 15.0;

		for (ArrayList<Double> point : points) {

			for (int j = 0; j < point.size(); j++) {

				if (j % 2 == 0) // x
					point.set(j, point.get(j) + xShift);

				else // y
					point.set(j, point.get(j) + yShift);
			}
		}

		width = xMax + Math.abs(xMin) + 15.0;
		height = yMax + Math.abs(yMin) + 15.0;
		this.resize(width, height);

		for (int i = 0; i < centersCoords.length; i++) {
			hexagons[i] = new Hexagon2(hexagonsCoords[i], points.get(i));
		}
	}

	private ArrayList<Double> getHexagonPoints(double xCenter, double yCenter) {

		ArrayList<Double> points = new ArrayList<>();

		points.add(xCenter);
		points.add(yCenter - 29.5);

		points.add(xCenter + 26.0);
		points.add(yCenter - 14.0);

		points.add(xCenter + 26.0);
		points.add(yCenter + 14.0);

		points.add(xCenter);
		points.add(yCenter + 29.5);

		points.add(xCenter - 26.0);
		points.add(yCenter + 14.0);

		points.add(xCenter - 26.0);
		points.add(yCenter - 14.0);

		return points;
	}

	public double getWidth() {
		return width;
	}

	public double getHeight() {
		return height;
	}

	protected Text createText(String string, double fontSize) {
		Text text = new Text(string);
		text.setBoundsType(TextBoundsType.VISUAL);
		text.setStyle("-fx-font-family: \"Verdana\";" + "-fx-font-size: " + fontSize + "px;");
		return text;
	}

	protected Text createText(String string) { //Pour compatibilité si l'ancienne signature est utilisée ailleurs
		return createText(string, 10); // Taille par défaut pour les IDs
	}

	protected void drawHexagons() {

		// Effacer les anciens textes avant de redessiner
		this.getChildren().removeAll(texts);
		this.getChildren().removeAll(cycleTypeTexts);
		texts.clear();
		cycleTypeTexts.clear();
		System.out.println("MoleculeGroup.drawHexagons() - Molecule: " + (molecule != null ? molecule.toString() : "null") + ", Nb Hexagons: " + hexagons.length);


		for (int i = 0; i < hexagons.length; i++) {
			Hexagon2 hexagon = hexagons[i];
			this.getChildren().add(hexagon);

			if (writeText) {
				// Texte pour l'ID de l'hexagone
				Text idText = createText(Integer.toString(i), 10); // Taille de police réduite

				// Positionner le texte de l'ID (par exemple, en haut au centre de l'hexagone)
				double idTextX = centersCoords[i].getX() + xShift - idText.getLayoutBounds().getWidth() / 2;
				double idTextY = centersCoords[i].getY() + yShift - 3; // Ajuster pour le positionnement
				idText.setX(idTextX);
				idText.setY(idTextY);
				texts.add(idText);
				this.getChildren().add(idText);

				// Texte pour le type de cycle de l'hexagone
				CycleType cycleType = molecule.getHexagonCycleType(i);
				if (cycleType != CycleType.UNKNOWN && !cycleType.toString().isEmpty()) { // N'afficher que si le type est connu et pertinent
					Text cycleTypeText = createText(cycleType.toString(), 9); // Taille de police encore plus réduite
					// Positionner le texte du type de cycle (par exemple, en bas au centre de l'hexagone)
					double cycleTypeTextX = centersCoords[i].getX() + xShift - cycleTypeText.getLayoutBounds().getWidth() / 2;
					double cycleTypeTextY = centersCoords[i].getY() + yShift + 10; // Ajuster pour le positionnement
					cycleTypeText.setX(cycleTypeTextX);
					cycleTypeText.setY(cycleTypeTextY);
					cycleTypeTexts.add(cycleTypeText);
					this.getChildren().add(cycleTypeText);
				}
			}
		}
	}

	public void removeTexts() {
		this.getChildren().removeAll(texts);
		texts.clear();
		if (cycleTypeTexts != null) {
			this.getChildren().removeAll(cycleTypeTexts);
			cycleTypeTexts.clear();
		}
	}

}