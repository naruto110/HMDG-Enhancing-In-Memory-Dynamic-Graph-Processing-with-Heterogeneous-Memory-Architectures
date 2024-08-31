#include<iostream>
using namespace std;
 
 
typedef size_t ElenmentType; 
//平衡二叉树的结构 
typedef struct AVLNode{
	int depth;//深度
    int tombstone; // 0:delete 1:insert
	struct AVLNode *left;//左孩子 
	struct AVLNode *right;//右孩子 
	struct AVLNode *parent;//父结点 
	ElenmentType value; //值 
}AVLtree,Tree;
//初始化 
Tree* avlInit(ElenmentType value){
	Tree* root=new Tree();
	root->parent = NULL;
    root->depth = 0;
    root->left = root->right = NULL;
    root->value=value; 
    root->tombstone = 1; // insertion
    return root;
}

Tree* avlInit(ElenmentType value, int tombstone){
	Tree* root=new Tree();
	root->parent = NULL;
    root->depth = 0;
    root->left = root->right = NULL;
    root->value=value; 
    root->tombstone = tombstone; 
    return root;
}
 
//LL型调整函数,执行右旋
Tree* LL_rotate(Tree *root){
    Tree *temp;
	temp = root->left;
	root->left = temp->right;
	temp->right = root;
	return temp; 
   
}
//RR型调整函数,执行左旋 
Tree* RR_rotate(Tree * root){
    Tree* temp;
    temp = root->right;
    root->right = temp->left;
    temp->left = root;
    return temp;
}
//LR型调整函数，先左旋转，再右旋转
Tree* LR_rotate(Tree* root){
	Tree* temp;
	temp = root->left;
    root->left =RR_rotate(temp);
    return LL_rotate(root);
} 
//RL型调整函数，先右旋转，再左旋转
Tree* RL_rotate(Tree* root){
	 Tree* temp;
	 temp = root->right;
    root->right=LL_rotate(temp);
    return RR_rotate(root);
} 
 
//树高
int height(Tree* root)
{
    if (root == NULL)
        return 0;
    int max;
    int left=height(root->left);
    int right=height(root->right);
    if(left>=right)
    max=left;
    else
    max=right;
    return max+1;
}
//求叶子节点个数
int  GetSumOfLeafNode(Tree* root)
{
	if(root == NULL)
		return 0;
	
	if(root->left == NULL && root->right == NULL)
		return 1;
	else
	{
		return GetSumOfLeafNode(root->left) 
			+ GetSumOfLeafNode(root->right);
	}
}
//平衡因子，即当前节点左右子树的差
int diff(Tree* root)
{
    return height(root->left) - height(root->right);
}
 
//平衡操作
Tree* avlBalance(Tree* root)
{
    int balanceFactor = diff(root);//平衡因子（左右子树高度差）
    if (balanceFactor > 1)//左子树高于右子树
    {
        if (diff(root->left) > 0)
		//LL的情况 
            root=LL_rotate(root);
        else
		//LR的情况 
            root=LR_rotate(root);
    }
    else if (balanceFactor < -1)//右子树高于左子树
    {
        if (diff(root->right) > 0)
		//RL的情况 
            root=RL_rotate(root);
        else
		//RR的情况 
            root=RR_rotate(root);
    }
    return root;
} 
//插入结点
Tree* avlInsert(Tree* root, ElenmentType k, int tombstone)
{
    if (NULL == root)
    {
        root = avlInit(k, tombstone);//如果根结点为null，则直接将值为根结点 
        if(root==NULL)
        	cout<<"creat failed"<<endl; 
        return root;
    }
    else if (k < root->value)
    {
        root->left = avlInsert(root->left, k, tombstone);//递归左子树
        root = avlBalance(root);//平衡操作
    }
    else if (k>root->value)
    {
        root->right = avlInsert(root->right, k, tombstone);//递归右子树
        root = avlBalance(root);//平衡操作
    }
    return root;
} 
//前序遍历
void displayTree_Pre(Tree* node){
        if(node == NULL) return;
        cout<<node->value<<" ";
        if(node->left != NULL){
            displayTree_Pre(node->left);
        }
        if(node->right != NULL){
            displayTree_Pre(node->right);
        }
}
//中序遍历
void displayTree_Mid(Tree* node){
        if(node == NULL) return;
        if(node->left != NULL){
            displayTree_Mid(node->left);
        }
        cout<<node->value<<" ";
        if(node->right != NULL){
            displayTree_Mid(node->right);
        }
}

//中序遍历
// src的作用是隔离非连续区域（假设不存在自环）
void displayTree_Mid(Tree* node, vector<size_t> &vec, size_t src){
    if(node == NULL) return;    
    if(node->left != NULL){
        displayTree_Mid(node->left, vec, src);
    }

    // if(src == 9858)
    //     cout<<node->value<<" ";

    if (vec.size() >= 1) // 不止一个元素
    {
        if ( (vec[vec.size()-1] != src) && (node->value - vec[vec.size()-1] > 1) )
        {
            vec.push_back(src);  // 当前待插入值相对最后一个元素值差值大于1则不连续，插入src用以标记
        }
    }
    vec.push_back(node->value);
    if(node->right != NULL){
        displayTree_Mid(node->right, vec, src);
    }
}

//中序遍历
// src的作用是隔离非连续区域（假设不存在自环）
void displayTree_Mid(Tree* node, size_t src, size_t &edge_num){
        if(node == NULL) return;
        if(node->left != NULL){
            displayTree_Mid(node->left, src, edge_num);
        }
        // cout<<node->value<<" ";
        // if (vec.size() >= 1) // 不止一个元素
        // {
        //     if ( (vec[vec.size()-1] != src) && (node->value - vec[vec.size()-1] > 1) )
        //     {
        //         vec.push_back(src);  // 当前待插入值相对最后一个元素值差值大于1则不连续，插入src用以标记
        //     }
        // }
        // vec.push_back(node->value);
        edge_num++;

        if(node->right != NULL){
            displayTree_Mid(node->right, src, edge_num);
        }
}

//后序遍历
void displayTree_Post(Tree* node){
        if(node == NULL) return;
        if(node->left != NULL){
            displayTree_Post(node->left);
        }
        if(node->right != NULL){
            displayTree_Post(node->right);
        }
        cout<<node->value<<" ";
}
//查找value 
Tree* binaryTreeSearch(Tree *node,int value){
	if(node->value==value)
		return node;
	//大于，在左边找 
	else if(node->value>value){
		if(node->left!=NULL)
			return binaryTreeSearch(node->left,value);
		else return NULL;
	}
	//否则，在右边找 
	else{
		if(node->right!=NULL)
			return binaryTreeSearch(node->right,value);
		else
			return NULL;
	}
}
 
//平衡二叉树最大值 
ElenmentType tree_max(Tree *node){
	int value;
	value=node->value;
	if(node->right!=NULL)
	return tree_max(node->right);
	else
	return value;
	
}
//平衡二叉树最小值 
ElenmentType tree_min(Tree *node){
	int value;
	value=node->value;
	if(node->left!=NULL)
	return tree_min(node->left);
	else
	return value;
}
 
//删除结点
Tree* avlDelete(Tree *root, const ElenmentType k)
{
    if (NULL == root)
        return root;
    if (!binaryTreeSearch(root,k))//查找删除元素是否存在
    {
        cout<<"Delete failed"<<endl;
        return root;
    }
 
    if (k == root->value)//根节点
    {
        if (root->left!=NULL&&root->right!=NULL)//左右子树都非空
        {
            if (diff(root) > 0)//左子树更高，在左边删除
            {
                root->value = tree_max(root->left);//以左子树的最大值替换当前值
                root->left = avlDelete(root->left, root->value);//删除左子树中已经替换上去的节点
            }
            else//右子树更高，在右边删除
            {
                root->value = tree_min(root->right);//以右子树的最小值替换当前值
                root->right = avlDelete(root->right, root->value);//删除右子树中已经替换上去的节点
            }
        }
        else//有一个孩子、叶子节点的情况合并
        {
                Tree * tmp = root;
                root = (root->left) ? (root->left) :( root->right);
                delete tmp;
                tmp = NULL;
        }
    }
    //往左边删 
    else if (k < root->value)
    {
        root->left = avlDelete(root->left, k);
        //不满足平衡条件
		if (diff(root) < -1)
        {
            if (diff(root->right) > 0)
            {
                root = RL_rotate(root);
            }
            else
            {
                root = RR_rotate(root);
            }
        }
    }
    //往右边删 
    else
    {
        root->right = avlDelete(root->right, k);
        //不满足平衡 条件 
        if (diff(root) > 1)
        {
            if (diff(root->left) < 0)
            {
                root = LR_rotate(root);
            }
            else
            {
                root = LL_rotate(root);
            }
        }
    }
    return root;
}
 
 
 
// int main(){
// 	int a[10]={3,5,4,2,6,8,10,1,7,9};
// 	Tree *root=NULL;
// 	for(int i=0;i<10;i++){
// 		root = Insert(root,a[i]);
// 	}
// 	cout<<"Preorder Traversal: ";
// 	displayTree_Pre(root);
// 	cout<<endl; 
// 	cout<<"Midorder Traversal: ";
// 	displayTree_Mid(root);
// 	cout<<endl; 
// 	cout<<"Tree Height:"<<height(root)<<endl; 
// 	int value = 33;
// 	cout<<"Max key: "<<tree_max(root)<<endl;
// 	cout<<"Min key: "<<tree_min(root)<<endl;
// 	cout<<"Num of leaf: "<<GetSumOfLeafNode(root)<<endl; 
// 	Tree* obj;
// 	if((obj=binaryTreeSearch(root,value))==NULL){
// 		cout<<value<<"No"<<endl;
// 	}
// 	else 
// 	    cout<<value<<"Yes"<<endl;
// 	root = Insert(root,5);
// 	cout<<"insert 5, Preorder Traversal: ";
//     displayTree_Pre(root);
// 	cout<<endl;
// 	cout<<"insert 5, Midorder Traversal: ";
// 	displayTree_Mid(root);
// 	cout<<endl;
// 	int del=79;
// 	if (!binaryTreeSearch(root,del))//查找删除元素是否存在
//     {
//         cout<<"Delete failed, no result"<<endl;
//     }else{
//     	Delete(root,del);
//     	cout<<"Delete node, Preorder Traversal: ";
//     	displayTree_Pre(root);
//     	cout<<endl;
//     	cout<<"Delete node, Midorder Traversal: ";
//     	displayTree_Mid(root);
//     	cout<<endl;
// 	} 
// 	return 0;
// }