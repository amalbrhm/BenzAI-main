package generator.properties.model;

import constraints.PentagonNumberConstraint;
import generator.properties.model.filters.PentagonNumberFilter;
import view.generator.ChoiceBoxCriterion;
import view.generator.boxes.HBoxModelCriterion;
import view.generator.boxes.HBoxPentagonNumberCriterion;
import view.primaryStage.ScrollPaneWithPropertyList;

public class PentagonNumberProperty extends ModelProperty {

    public PentagonNumberProperty() {
        super("pentagons", "Number of pentagons", new PentagonNumberConstraint(), new PentagonNumberFilter());
    }

    @Override
    public HBoxModelCriterion makeHBoxCriterion(ScrollPaneWithPropertyList parent, ChoiceBoxCriterion choiceBoxCriterion) {
        return new HBoxPentagonNumberCriterion(parent, choiceBoxCriterion);
    }

    @Override
    public int computeHexagonNumberUpperBound() {
        return Integer.MAX_VALUE;
    }

    @Override
    public int computeNbCrowns() {
        return ((ModelPropertySet) this.getPropertySet()).getHexagonNumberUpperBound() > 0 ?
                (((ModelPropertySet) this.getPropertySet()).getHexagonNumberUpperBound() + 2) / 2 :
                Integer.MAX_VALUE;
    }
}