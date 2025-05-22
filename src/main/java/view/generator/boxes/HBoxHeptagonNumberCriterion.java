package view.generator.boxes;

import generator.properties.model.ModelPropertySet;
import generator.properties.model.expression.BinaryNumericalExpression;
import view.generator.ChoiceBoxCriterion;
import view.primaryStage.ScrollPaneWithPropertyList;

public class HBoxHeptagonNumberCriterion extends HBoxBoundingCriterion {

    public HBoxHeptagonNumberCriterion(ScrollPaneWithPropertyList parent, ChoiceBoxCriterion choiceBoxCriterion) {
        super(parent, choiceBoxCriterion);
    }

    @Override
    public void addPropertyExpression(ModelPropertySet modelPropertySet) {
        if (isValid())
            modelPropertySet.getById("heptagons").addExpression(new BinaryNumericalExpression("heptagons", getOperatorChoiceBox().getValue(), Integer.decode(getFieldValue().getText())));
    }
}