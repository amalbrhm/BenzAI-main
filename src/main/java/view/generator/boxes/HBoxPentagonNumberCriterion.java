package view.generator.boxes;

import generator.properties.model.ModelPropertySet;
import generator.properties.model.expression.BinaryNumericalExpression;
import view.generator.ChoiceBoxCriterion;
import view.primaryStage.ScrollPaneWithPropertyList;

public class HBoxPentagonNumberCriterion extends HBoxBoundingCriterion {

    public HBoxPentagonNumberCriterion(ScrollPaneWithPropertyList parent, ChoiceBoxCriterion choiceBoxCriterion) {
        super(parent, choiceBoxCriterion);
    }

    @Override
    public void addPropertyExpression(ModelPropertySet modelPropertySet) {
        if (isValid())
            modelPropertySet.getById("pentagons").addExpression(new BinaryNumericalExpression("pentagons", getOperatorChoiceBox().getValue(), Integer.decode(getFieldValue().getText())));
    }
}