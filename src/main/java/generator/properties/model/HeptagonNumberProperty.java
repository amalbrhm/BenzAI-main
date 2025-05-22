package generator.properties.model;

import constraints.HeptagonNumberConstraint;
import generator.properties.model.filters.HeptagonNumberFilter;
import view.generator.ChoiceBoxCriterion;
import view.generator.boxes.HBoxModelCriterion;
import view.generator.boxes.HBoxHeptagonNumberCriterion;
import view.primaryStage.ScrollPaneWithPropertyList;

public class HeptagonNumberProperty extends ModelProperty {

    public HeptagonNumberProperty() {
        super("heptagons", "Number of heptagons", new HeptagonNumberConstraint(), new HeptagonNumberFilter());
    }

    @Override
    public HBoxModelCriterion makeHBoxCriterion(ScrollPaneWithPropertyList parent, ChoiceBoxCriterion choiceBoxCriterion) {
        return new HBoxHeptagonNumberCriterion(parent, choiceBoxCriterion);
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